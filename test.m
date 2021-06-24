clear;
clc;
profile on;
PathSet();
[meshType, nV, nv, numDecompose, fileName, epsArap, epsSchur] = SetParameter();
[Vertex, Face] = MeshGeneration(meshType, nV, fileName);
[nF, I, nI, B, nB, MC] = MeshInfo(Vertex, Face, nV);
[DS, DE, DW, Map, sepEdge] = Decomp(MC, Vertex, nV, numDecompose);

%%	solve Hessian system

%   solve Hessian system directly
[Hx, Hy, t1, G, dz1, Vertex0] = testfunc1(Vertex, Vertex, Face, 1);

%   compute the Schur matrix blocks
[dsInd, dwInd, deInd, dsIndall, dweInd, ide, dssize] = SchurSystemIndex(DS, DE, DW, numDecompose);
[MSSx, MSEx, MSEallx, MWWx, MWEx, MEEx] ...
= SchurSystemM(Hx, dsInd, dwInd, deInd, dsIndall, numDecompose);
[MSSy, MSEy, MSEally, MWWy, MWEy, MEEy] ...
= SchurSystemM(Hy, dsInd, dwInd, deInd, dsIndall, numDecompose);

%   compute the faces of wirebasket/edge
vwe = Vertex0(dweInd, :);
vweorigin = Vertex(dweInd, :);
[FinWE, finweInd] = FindFinWE(Face, dweInd, nV);
fwe = FindFforWE(FinWE, dweInd);

%   solve Hessian system with PCG
[iter, t2, dz2] = testfunc2(MSSx, MSEx, MSEallx, MWWx, MWEx, MEEx, MSSy, MSEy, MSEally, MWWy, MWEy, MEEy, ...
                            dsInd, dwInd, deInd, dssize, vwe, vweorigin, fwe, ide, G, ...
                            nV, numDecompose, epsSchur);
                        
%   the error
error = RMSError1(dz1, dz2);

%%
profile viewer;