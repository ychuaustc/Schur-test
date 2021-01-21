%   parameterization using local-global method (ARAP energy) and preconditioned conjugate solver


clear;
clc;


%   profile on;
addpath(genpath('NewtonParam'));
fprintf('Mesh Parameterization Test.\n\n\n');


%%  set parameters
[meshType, nV, numDecomposeTemp, nv, wirebasketType, solverType, epsArap, epsSchur] = SetParameter();


%%  mesh generation
fprintf('Part I: Mesh generation\n\n');
[Vertex, Face] = MeshGeneration(meshType, nV);


%%  mesh information
[nF, I, nI, B, nB, MC] = MeshInfo(Vertex, Face, nV);


%%	mesh decomposition
fprintf('Part II: Mesh decomposition\n\n');
[DV, DS, DB, DE, DW, DWB, numDS, numDecompose] = Decompose(nV, nv, B, MC, wirebasketType, Vertex, Face, nF);


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

%   set the constraint
fixedP = [];    % fixed points
vecWOfixedP = setdiff(1:nV, fixedP)';
MWOF = M(vecWOfixedP, vecWOfixedP); % the coefficient matrix without the fixed points

%   compute the Schur matrix blocks
[MSS, MSE, MWW, MWE, MEE, dsInd, dwInd, deInd] = SchurSystemM(M, DS, DE, DW, numDecompose, fixedP);

%   compute the Larangian multiplier for the singular system
[LW, LE] = LagMultip(FXV, FYV, DW, DWB, dwInd, deInd, Face, nV, nF);

%   compute the wirebasket matrix block for the preconditioner
if solverType == 3
    MEE_W = MWirebasket(FCotTheta, Face, DW, DWB, deInd, nV, nF);
end



%%	Parameterization process
fprintf('Part IV: Parameterization using the local-global method\n\n');

%   compute the rotation L for the local phase (the first iteration step)
[FLVX, FLVY, MThetaLVX, MThetaLVY] = ArapL(FXV, FYV, FXU, FYU, FCotTheta, Face, nV, nF);

%   compute the ARAP energy E(U, L) for the global phase (the first iteration step)
EnUL = EnergyUL(FCotTheta, FLVX, FLVY, FXU, FYU);
EnULOld = EnUL;
EnULNew = 0;

%   the iteration
iterARAP = 0; % iteration count

while 1
	iterARAP = iterARAP + 1;
	CX = ArapSystemC(MThetaLVX);
	CY = ArapSystemC(MThetaLVY);
    
    switch solverType
        case 1  % direct solver
            CXWOF = CX(vecWOfixedP) - M(vecWOfixedP, fixedP) * VertexU(fixedP, 1);
            CYWOF = CY(vecWOfixedP) - M(vecWOfixedP, fixedP) * VertexU(fixedP, 2);
            [XUWOF, t_X{iterARAP}] = DirectSolver(MWOF, CXWOF);
            [YUWOF, t_Y{iterARAP}] = DirectSolver(MWOF, CYWOF);
            XU(vecWOfixedP) = XUWOF;
            YU(vecWOfixedP) = YUWOF;
        case 2  % conjugate gradient solver
            CX(vecWOfixedP) = CX(vecWOfixedP) - M(vecWOfixedP, fixedP) * VertexU(fixedP, 1);
            CY(vecWOfixedP) = CY(vecWOfixedP) - M(vecWOfixedP, fixedP) * VertexU(fixedP, 2);
            [CSX, CWX, bX] = SchurSystemC(MSS, MSE, MWW, MWE, CX, dsInd, dwInd, deInd, numDecompose);
            [CSY, CWY, bY] = SchurSystemC(MSS, MSE, MWW, MWE, CY, dsInd, dwInd, deInd, numDecompose);
            [XUWOF, iter_X{iterARAP}, t_X{iterARAP}] ...
            = SchurConjSolver(MSS, MSE, MWW, MWE, MEE, CSX, CWX, bX, dsInd, dwInd, deInd, nV - 1, numDecompose, epsSchur);
            [YUWOF, iter_Y{iterARAP}, t_Y{iterARAP}]...
            = SchurConjSolver(MSS, MSE, MWW, MWE, MEE, CSY, CWY, bY, dsInd, dwInd, deInd, nV - 1, numDecompose, epsSchur);
            XU(vecWOfixedP) = XUWOF(vecWOfixedP);
            YU(vecWOfixedP) = YUWOF(vecWOfixedP);
        case 3  % preconditioned conjugate gradient solver
            CX(vecWOfixedP) = CX(vecWOfixedP) - M(vecWOfixedP, fixedP) * VertexU(fixedP, 1);
            CY(vecWOfixedP) = CY(vecWOfixedP) - M(vecWOfixedP, fixedP) * VertexU(fixedP, 2);
            [CSX, CWX, bX] = SchurSystemC(MSS, MSE, MWW, MWE, CX, dsInd, dwInd, deInd, numDecompose);
            [CSY, CWY, bY] = SchurSystemC(MSS, MSE, MWW, MWE, CY, dsInd, dwInd, deInd, numDecompose);
            [XUWOF, iter_X{iterARAP}, t_X{iterARAP}] ...
            = SchurConjPreSolver(MSS, MSE, MWW, MWE, MEE, MEE_W, LW, LE, CSX, CWX, bX, dsInd, dwInd, deInd, nV - 1, numDecompose, epsSchur);
            [YUWOF, iter_Y{iterARAP}, t_Y{iterARAP}]...
            = SchurConjPreSolver(MSS, MSE, MWW, MWE, MEE, MEE_W, LW, LE, CSY, CWY, bY, dsInd, dwInd, deInd, nV - 1, numDecompose, epsSchur);
            XU(vecWOfixedP) = XUWOF(vecWOfixedP);
            YU(vecWOfixedP) = YUWOF(vecWOfixedP);
        otherwise
            quit(1)
    end
    
	VertexU = [XU, YU];
	[FXU, FYU] = FaceCoord(VertexU, Face);
	[FLVX, FLVY, MThetaLVX, MThetaLVY] = ArapL(FXV, FYV, FXU, FYU, FCotTheta, Face, nV, nF);    % compute the rotation L for
                                                                                                % current iteration
 	EnULNew = EnergyUL(FCotTheta, FLVX, FLVY, FXU, FYU);    % compute the ARAP energy for current iteration
    
 	if abs(EnULNew - EnULOld) / nF >= epsArap
        EnULOld = EnULNew;
    else
        break;
	end
end

%   display the result
if solverType == 2 || solverType == 3
    iterSchur = 0;
    for i = 1:iterARAP
        iterSchur = iterSchur + iter_X{i} + iter_Y{i};
    end
    iterSchur = floor(iterSchur / (iterARAP * 2));
    fprintf('local-global iteration steps: %d \n', iterARAP);
    fprintf('PCG iteration steps (in each global phase): %d \n', iterSchur);
end
VertexU = full(VertexU);
DrawMesh(VertexU, Face, nV);


% profile viewer;