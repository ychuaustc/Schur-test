%%  check the validation of all solvers by comparing them to direct solver with given M and CT (tuple)
%%
%%	Input:
%%          nV:         mesh size
%%          nv:         size of each decomposition
%%          epsSchur:   value for Schur solver convergence check
%%
%%	Output:

function [] = SolverCheck(nV, nv, epsSchur)

%%  construction of M and C (using posisson equation and p1 element on structured mesh)
%   mesh generation
[Vertex, Face] = GenerateST(nV);
%   mesh information and decomposition
[nF, I, nI, B, nB, MC] = MeshInfoUT(Vertex, Face, nV);
[DV, DS, DB, DE, DW, DWB, numDS, numDecompose] = DecomposeUT(nV, nv, B, MC);
%   mesh display
% DrawMesh(Vertex, Face, nV);
% DrawDecomposedMesh(Vertex, Face, DV, nV, numDecompose);
%   M construction
M = Poisson1D_P1_M(I, nI, B, nV);
%   C construction
nC = 2; %   size of tuple CT
for i = 1:nC
    CT{1, i} = rand(nV, 1);
end

%%  solving MX = C with different solvers
%   constructing Schur system
[MSS, MSE, MWW, MWE, MEE, CST, CWT, S, bT, dsInd, dwInd, deInd] = SchurSystem(M, CT, DS, DB, DE, DW, DWB, numDecompose);
%   direct solver
for i = 1:nC
    [XT_D{i}, t_D{i}] = DirectSolver(M, CT{i});
end
%   Schur solver
for i = 1:nC
    [XT_S{i}, t_S{i}] = SchurSolver(MSS, MSE, MWW, MWE, CST{i}, CWT{i}, S, bT{i}, dsInd, dwInd, deInd, nV, numDecompose);
end
%   Schur conjugate solver
for i = 1:nC
    [XT_SC{i}, iterT_SC{i}, t_SC{i}] = SchurConjSolver(MSS, MSE, MWW, MWE, MEE, CST{i}, CWT{i}, bT{i}, dsInd, dwInd, deInd, nV, numDecompose, epsSchur);
end
%   Schur conjugate solver using preconditioner
for i = 1:nC
    [XT_SCP{i}, iterT_SCP{i}, t_SCP{i}] = SchurConjPreSolver(MSS, MSE, MWW, MWE, MEE, CST{i}, CWT{i}, bT{i}, dsInd, dwInd, deInd, nV, numDecompose, epsSchur);
end

%%  solver comparision using L1-error
%   Schur solver
err_S = 0;
for i = 1:nC
    err_S = err_S + L1error(XT_D{i}, XT_S{i}, nV);
end
err_S = err_S / nC;
%   Schur conjugate solver
err_SC = 0;
for i = 1:nC
    err_SC = err_SC + L1error(XT_D{i}, XT_SC{i}, nV);
end
err_SC = err_SC / nC;
%   Schur conjugate solver using preconditioner
err_SCP = 0;
for i = 1:nC
    err_SCP = err_SCP + L1error(XT_D{i}, XT_SCP{i}, nV);
end
err_SCP = err_SCP / nC;

%%
end