%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Yuchen Hua, USTC, hyc12908@ustc.edu.cn
%%          
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear;
clc;
%%
profile on;
%%
addpath(genpath('D:\OneDrive\USTC\USTC\research\projects\SchurParallel\SchurBasedSolver - new\NewtonParam'));
%%
fprintf('Mesh Parameterization.\n\n\n\n');
%%
epsArap = 1e-7;     %   value for ARAP solver convergence check
epsSchur = 1e-6;	%   value for Schur solver convergence check
%%
solver = 3;         %   1: direct solver; 2: conjugate gradient solver; 3: preconditioned conjugate gradient solver
meshType = 1;       %   1: structured triangular mesh; 2: unstructured triangular mesh; 3: read meshes from a "obj" file
wirebsketType = 1;  %   1: based on intersection points; 2: based on edges
m1 = 3;             %   parameter for mesh size
m2 = 1;             %   parameter for decomposition number
%%
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  mesh treatment
%%
%%  mesh generation / read
%   mesh generation
fprintf('Part I: Mesh generation\n\n');
nV = (2^m1 + 1)^2;              %	mesh size
numDecomposeTemp = 2^(2 * m2);  %	expected number of decompositions
nv = (2^(m1 - m2) + 1)^2;       %	expected size of each decomposition
if meshType == 1
    [Vertex, Face] = GenerateST(nV);
end
if meshType == 2
    [Vertex, Face] = GenerateUT(nV);
end
if meshType == 3
    [Vertex, Face] = ...
    readObj('NewtonParam\camelhead_slim');
end
VertexX = Vertex(:, 1);
VertexY = Vertex(:, 2);
%   mesh read
%%  mesh information
[nF, I, nI, B, nB, MC] = MeshInfoUT(Vertex, Face, nV);
%% mesh decomposition
fprintf('Part II: Mesh decomposition\n\n');
% [DV, ~, ~, ~, ~, ~, numDS, numDecompose] = Decompose(nV, nv, B, MC);
% [DS, DB, DE, DW, DWB] = ModifyDecompose(DV, numDecompose, MC, nV);
[DV, DS, DE, DW, DWB, numDecompose] = DecomposeTest(nV);
%%  restriction
% dsInd = find(DS{1})';
% sqnDV = sqrt(size(find(DV{1}), 1));
% fixedP = dsInd((sqnDV - 1) * (sqnDV - 1) / 2 + (sqnDV + 1) / 2)';
% fixedP = dsInd(10)';
fixedP = [];
vecWOfixedP = setdiff(1:nV, fixedP)';
%%
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  parameteization using arap energy
%   geometry information of original mesh
[FXV, FYV] = FaceCoord(Vertex, Face);
[FCotTheta, MCotTheta] = GeoMeshInfo(Face, FXV, FYV, nV);
%	Lagrangian multiplier
[LW, LE] = LagMultip(FXV, FYV, DW, DWB, DE, Face, nV, nF);
%   initial parameterization
VertexU = Initialization(MC, I, B, nB, nV);
% DrawMesh(VertexU, Face, nV);
[FXU, FYU] = FaceCoord(VertexU, Face);
%   compute parameterization using arap method
%   compute the rotation L
[FLVX, FLVY, MThetaLVX, MThetaLVY] = ArapL(FXV, FYV, FXU, FYU, FCotTheta, Face, nV, nF);
%   compute the energy E(U, L)
EnUL = EnergyUL(FCotTheta, FLVX, FLVY, FXU, FYU);
EnULOld = EnUL;
EnULNew = 0;
%%	iteration steps
%	compute M for ARAP system MX = C
M = ArapSystemM(MCotTheta, nV);
%% deal with restriction
%   M
MWOF = M(vecWOfixedP, vecWOfixedP);
%
XU = VertexU(:, 1);
YU = VertexU(:, 2);
%   decomposition
DS{1}(fixedP) = 0;
[MSS, MSE, MWW, MWE, MEE, dsInd, dwInd, deInd] = SchurSystemM(M, DS, DE, DW, numDecompose);
%
iterARAP_D = 0;
%   direct solver
if solver == 1
    while 1
        iterARAP_D = iterARAP_D + 1;
        CX = ArapSystemC(MThetaLVX);
        CY = ArapSystemC(MThetaLVY);
        CXWOF = CX(vecWOfixedP) - M(vecWOfixedP, fixedP) * VertexU(fixedP, 1);
        CYWOF = CY(vecWOfixedP) - M(vecWOfixedP, fixedP) * VertexU(fixedP, 2);
        [XUWOF, t_DX{iterARAP_D}] = DirectSolver(MWOF, CXWOF);
        [YUWOF, t_DY{iterARAP_D}] = DirectSolver(MWOF, CYWOF);
        XU(vecWOfixedP) = XUWOF;
        YU(vecWOfixedP) = YUWOF;
        VertexU = [XU, YU];
        [FXU, FYU] = FaceCoord(VertexU, Face);
        [FLVX, FLVY, MThetaLVX, MThetaLVY] = ArapL(FXV, FYV, FXU, FYU, FCotTheta, Face, nV, nF);
        EnULNew = EnergyUL(FCotTheta, FLVX, FLVY, FXU, FYU);
        if abs(EnULNew - EnULOld) / nF >= epsArap
            EnULOld = EnULNew;
        else
            break;
        end
    end
end
%   Schur conjugate solver
if solver == 2
    while 1
        iterARAP_D = iterARAP_D + 1;
        CX = ArapSystemC(MThetaLVX);
        CY = ArapSystemC(MThetaLVY);
        CX(vecWOfixedP) = CX(vecWOfixedP) - M(vecWOfixedP, fixedP) * VertexU(fixedP, 1);
        CY(vecWOfixedP) = CY(vecWOfixedP) - M(vecWOfixedP, fixedP) * VertexU(fixedP, 2);
        [CSX, CWX, bX] = SchurSystemC(MSS, MSE, MWW, MWE, CX, dsInd, dwInd, deInd, numDecompose);
        [CSY, CWY, bY] = SchurSystemC(MSS, MSE, MWW, MWE, CY, dsInd, dwInd, deInd, numDecompose);
        [XUWOF, iterX_SC{iterARAP_D}, tX_SC{iterARAP_D}] ...
        = SchurConjSolver(MSS, MSE, MWW, MWE, MEE, CSX, CWX, bX, dsInd, dwInd, deInd, nV - 1, numDecompose, epsSchur);
        [YUWOF, iterY_SC{iterARAP_D}, tY_SC{iterARAP_D}]...
        = SchurConjSolver(MSS, MSE, MWW, MWE, MEE, CSY, CWY, bY, dsInd, dwInd, deInd, nV - 1, numDecompose, epsSchur);
        XU(vecWOfixedP) = XUWOF(vecWOfixedP);
        YU(vecWOfixedP) = YUWOF(vecWOfixedP);
        VertexU = [XU, YU];
        [FXU, FYU] = FaceCoord(VertexU, Face);
        [FLVX, FLVY, MThetaLVX, MThetaLVY] = ArapL(FXV, FYV, FXU, FYU, FCotTheta, Face, nV, nF);
        EnULNew = EnergyUL(FCotTheta, FLVX, FLVY, FXU, FYU);
        if abs(EnULNew - EnULOld) / nF >= epsArap
            EnULOld = EnULNew;
        else
            break;
        end
    end
end
%   Schur conjugate solver using preconditioner
if solver == 3
    MEE_W = MWirebasket(FCotTheta, Face, DW, DWB, deInd, nV, nF);
    while 1
        iterARAP_D = iterARAP_D + 1;
        CX = ArapSystemC(MThetaLVX);
        CY = ArapSystemC(MThetaLVY);
        CX(vecWOfixedP) = CX(vecWOfixedP) - M(vecWOfixedP, fixedP) * VertexU(fixedP, 1);
        CY(vecWOfixedP) = CY(vecWOfixedP) - M(vecWOfixedP, fixedP) * VertexU(fixedP, 2);
        [CSX, CWX, bX] = SchurSystemC(MSS, MSE, MWW, MWE, CX, dsInd, dwInd, deInd, numDecompose);
        [CSY, CWY, bY] = SchurSystemC(MSS, MSE, MWW, MWE, CY, dsInd, dwInd, deInd, numDecompose);
        [XUWOF, iterX_SCP{iterARAP_D}, tX_SCP{iterARAP_D}] ...
        = SchurConjPreSolver(MSS, MSE, MWW, MWE, MEE, MEE_W, LW, LE, CSX, CWX, bX, dsInd, dwInd, deInd, nV - 1, numDecompose, epsSchur);
        [YUWOF, iterY_SCP{iterARAP_D}, tY_SCP{iterARAP_D}]...
        = SchurConjPreSolver(MSS, MSE, MWW, MWE, MEE, MEE_W, LW, LE, CSY, CWY, bY, dsInd, dwInd, deInd, nV - 1, numDecompose, epsSchur);
        XU(vecWOfixedP) = XUWOF(vecWOfixedP);
        YU(vecWOfixedP) = YUWOF(vecWOfixedP);
        VertexU = [XU, YU];
        [FXU, FYU] = FaceCoord(VertexU, Face);
        [FLVX, FLVY, MThetaLVX, MThetaLVY] = ArapL(FXV, FYV, FXU, FYU, FCotTheta, Face, nV, nF);
        EnULNew = EnergyUL(FCotTheta, FLVX, FLVY, FXU, FYU);
        if abs(EnULNew - EnULOld) / nF >= epsArap
            EnULOld = EnULNew;
        else
            break;
        end
    end
end
%
VertexU = full(VertexU);
%%
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  show result
%%  draw original mesh
DrawMesh(VertexU, Face, nV);
%%  draw decomposed mesh
% DrawDecomposedMesh(VertexU, Face, DV, nV, numDecompose);
%%
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% profile viewer;
%%
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%