clear;
clc;


addpath(genpath('NewtonParam'));
addpath(genpath('Hessian'));
addpath(genpath('H'));


[meshType, nV, numDecomposeTemp, nv, wirebasketType, solverType, epsArap, epsSchur, fileName] = SetParameter();
[Vertex, Face] = MeshGeneration(meshType, nV, fileName);
[nF, I, nI, B, nB, MC] = MeshInfo(Vertex, Face, nV);
[DV, DS, DB, DE, DW, DWB, numDS, numDecompose] = Decompose(nV, nv, B, MC, wirebasketType, Vertex, Face, nF);


M = getM(meshType, nV, fileName);


[MSS, MSE, MWW, MWE, MEE, dsInd, dwInd, deInd] = SchurSystemM(M, DS, DE, DW, numDecompose);


MEE_W = getH(meshType, nV, fileName, numDecompose);



XSol = rand(nV, 1);
C = M * XSol;
    
switch solverType
    case 1  % direct solver
        [XU, t] = DirectSolver(M, C);
    case 2  % conjugate gradient solver
        [CS, CW, b] = SchurSystemC(MSS, MSE, MWW, MWE, C, dsInd, dwInd, deInd, numDecompose);
        [XU, iter, t] ...
        = SchurConjSolver(MSS, MSE, MWW, MWE, MEE, CS, CW, b, dsInd, dwInd, deInd, nV, numDecompose, epsSchur);
    case 3  % preconditioned conjugate gradient solver
        [CS, CW, b] = SchurSystemC(MSS, MSE, MWW, MWE, C, dsInd, dwInd, deInd, numDecompose);
        [XU, iter, t] ...
        = SchurConjPreSolver(MSS, MSE, MWW, MWE, MEE, MEE_W, CS, CW, b, dsInd, dwInd, deInd, nV, numDecompose, epsSchur);
    otherwise
        quit(1)
end

if solverType == 2 || solverType == 3
    fprintf('The iteration time is %d\n', iter);
end

l2Error = norm(XSol - XU, 2) / nV;
fprintf('The L2-error of the result compared to the exact value is %d\n', l2Error);


% profile viewer;