%   parameterization using local-global method (ARAP energy) and preconditioned conjugate solver


clear;
clc;


profile on;
% delete(gcp('nocreate'));
% p = parpool;
PathSet();
fprintf('Mesh Parameterization Test.\n\n\n');


%%  set parameters
[meshType, nV, numDecompose, fileName, epsArap, epsSchur] = SetParameter();


%%  mesh generation
fprintf('Part I: Mesh generation\n\n');
[Vertex, Face] = MeshGeneration(meshType, nV, fileName);


%%  mesh information
[nF, I, nI, B, nB, MC] = MeshInfo(Vertex, Face, nV);


%%	mesh decomposition
fprintf('Part II: Mesh decomposition\n\n');
[DS, DE, DW] = Decomp(MC, Vertex, nV, numDecompose);
% [DS, DE, DW, numDecompose] = Decomp1(nV);


%%	parameterization pretreatment
fprintf('Part III: Pretreatment for the parameterization\n\n\n');

%   cotangent matrix of the mesh
[FXV, FYV] = FaceCoord(Vertex, Face);
[FCotTheta, MCotTheta] = CotMatrix(Face, FXV, FYV, nV);

%   initial parameterization
[VertexU, XU, YU] = Initialization(MC, I, B, nB, nV);
[FXU, FYU] = FaceCoord(VertexU, Face);

%	compute the coefficient matrix M for the ARAP system MX = C
M = ArapSystemM(MCotTheta, nV);

%   compute the Schur matrix blocks
[MSS, MSE, MWW, MWE, MEE, dsInd, dwInd, deInd] = SchurSystemM(M, DS, DE, DW, numDecompose);

%   compute the Larangian multiplier for the singular system
[LW, LE] = LagMultip(FXV, FYV, DW, DE, dwInd, deInd, Face, nV, nF);

%   compute the wirebasket matrix block for the preconditioner
MEE_W = MWirebasket(FCotTheta, Face, DW, DE, deInd, nV, nF);


%%	Parameterization process
fprintf('Part IV: Parameterization using the local-global method\n\n');

%   compute the rotation L for the local phase (the first iteration step)
[FLVX, FLVY, MThetaLVX, MThetaLVY] = ArapL(FXV, FYV, FXU, FYU, FCotTheta, Face, nV, nF);

%   compute the ARAP energy E(U, L) for the global phase (the first iteration step)
EnUL = EnergyUL(FCotTheta, FLVX, FLVY, FXU, FYU);
EnULOld = EnUL;

%   the iteration
%   direct solver
[VertexU_D, t_D, iterARAP_D] = IterDirect(MThetaLVX, MThetaLVY, M, Face, FXV, FYV, FCotTheta, EnULOld, nV, nF, epsArap);
%   conjugate gradient solver
[VertexU_C, t_C, iterARAP_C, iterSchur_C] = IterConj(MThetaLVX, MThetaLVY, MSS, MSE, MWW, MWE, MEE, ...
                                                     dsInd, dwInd, deInd, Face, FXV, FYV, FCotTheta, EnULOld, ...
                                                     nV, nF, numDecompose, epsArap, epsSchur);
%   preconditioned conjugate gradient solver
[VertexU_CP, t_CP, iterARAP_CP, iterSchur_CP] = IterConjPre(MThetaLVX, MThetaLVY, MSS, MSE, MWW, MWE, MEE, MEE_W, ...
                                                            LW, LE, dsInd, dwInd, deInd, Face, FXV, FYV, FCotTheta, ...
                                                            EnULOld, nV, nF, numDecompose, epsArap, epsSchur);


%%  display the result
ShowResult(iterARAP_D, t_D, iterARAP_C, iterSchur_C, t_C, iterARAP_CP, iterSchur_CP, t_CP)
for i = 2:4
    PlotMesh(VertexU_D, Face, i);
end


%%
% delete(gcp('nocreate'));
profile viewer;