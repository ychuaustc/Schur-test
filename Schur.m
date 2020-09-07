%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Yuchen Hua, USTC, hyc12908@ustc.edu.cn
%%
%%  parellel performance of the Schur compliment based algorithm on unstructured triangle meshes
%%
%%  The code is constructed with the following functions
%%      Main body:
%%          Schur.m:                main body
%%      Mesh treatment:
%%          GenerateS.m:            mesh generation (structured mesh)
%%          GenerateU.m:            mesh generation (unstructured mesh)
%%          MeshInfo.m:             information of mesh generated
%%      Mesh decomposition:
%%          Decompose.m:            return a decomposition of given mesh to Schur.m
%%          DecomposeByNeighbor.m:  return the update of vertice and connection matrix to Decompose.m
%%      Parallel computing test:
%%          Lsystem.m:              set up the linear system
%%          Parameterization.m:     parameterization by solving linear system
%%          ParameterizationD.m:    parameterization using conjugate solver
%%          ParameterizationC.m:	parameterization using conjugate solver
%%          ParameterizationCP.m:   parameterization using conjugate solver and preconditioner
%%      Display:
%%          DrawDecomposedMesh.m:   draw decomposed mesh
%%          DrawMeshCompare.m:      draw parameterized with original / decomposed mesh seperately
%%          
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear;
clc;

%%
profile on;

%%
fprintf('Schur compliment based solver for large scale linear systems.\n\n\n\n');

%% set up the parallel environment
mc = 2;
corenum = mc;                 % core number for parallel computing
% p = parpool(corenum);

%% mesh generation
fprintf('Part I: Mesh generation\n\n');
m1 = 7;
m2 = 3;
nv = (2^m1 + 1)^2;              % number of vertices
numDecomposeTemp = 2^(2 * m2);  % expected number of decompositions
n = (2^(m1 - m2) + 1)^2;        % expected number of vertice in each decomposition
% structured mesh
% [Vertex, Face] = GenerateS(nv);
% [Vertex, Face] = GenerateSR(nv);
% unstructured mesh
[Vertex, Face] = GenerateU(nv);
% mesh information
[nV, nF, B, nB, I, nI, MC] = MeshInfo(Vertex, Face);

%% mesh decomposition
fprintf('Part II: Mesh decomposition\n\n');
[DV, DB, ~, DS, DE, DW, DWB, numDecompose] = Decompose(n, nV, B, MC);

%% set up the linear system
fprintf('Part III: set up linear system\n\n');
[M, Cx, Cy, dsInd, dwInd, deInd, debInd, dewbInd, MS, MSE, MW, MWE, ME, MB, MSB, MWB, MWWB, CSx, CSy, CWx, CWy, CEx, CEy] ...
= Lsystem(Vertex, B, nB, I, nI, MC, DS, DE, DW, DB, DWB, numDecompose);

%% parameterization with original mesh
fprintf('Part IV: parameterization with original mesh\n\n');
[vertexo] = Parameterization(M, Cx, Cy, Vertex, I, B, nB);

%% parameterization with decomposed mesh
fprintf('Part V: parameterization with decomposed mesh\n\n');
% direct solver
[vertexdD] = ParameterizationD(dsInd, dwInd, deInd, MS, MSE, MW, MWE, ME, CSx, CSy, CWx, CWy, CEx, CEy, ...
                               Vertex, I, nI, B, nB, numDecompose);
% conjugate solver
[vertexdC] = ParameterizationC(dsInd, dwInd, deInd, MS, MSE, MW, MWE, ME, CSx, CSy, CWx, CWy, CEx, CEy, ...
                               Vertex, I, nI, B, nB, numDecompose);
% conjugate solver using preconditioner
[vertexdCP] = ParameterizationCP(dsInd, dwInd, deInd, debInd, dewbInd, MS, MSE, MW, MWE, ME, MB, MSB, MWB, MWWB, ...
                                 CSx, CSy, CWx, CWy, CEx, CEy, Vertex, I, nI, B, nB, numDecompose);

%% compute the l1 - error of vertex computed using different method
fprintf('Part IV: Error estimation: \n\n');
ed = l1error(vertexo, vertexdD, nV);
fprintf('L1-error for vertice obtained with original / decomposed mesh using direct solver: %f\n\n', ed);
ec = l1error(vertexo, vertexdC, nV);
fprintf('L1-error for vertice obtained with original / decomposed mesh using conjugate solver: %f\n\n', ec);
ecp = l1error(vertexo, vertexdCP, nV);
fprintf('L1-error for vertice obtained with original / decomposed mesh using conjugate solver with preconditioner: %f\n\n', ecp);

%% draw the result
% DrawDecomposedMesh(Vertex, Face, DV, nV, numDecompose);
DrawMeshCompare(vertexo, vertexdD, vertexdC, vertexdCP, Face);

%%
% delete(gcp('nocreate'));

%%
profile viewer;