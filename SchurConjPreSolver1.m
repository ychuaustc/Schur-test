function [X, iterT, t] = SchurConjPreSolver1(MSS, MSE, MSEall, MWW, MWE, MEE, MEE_W, LW, LE, CS, CW, b, dsInd, dwInd, deInd, MSSsolver, presolver, dssize, nV, numDecompose, epsSchur)
%
%   this function solves the Schur system using a preconditioned conjugate method
%
%   INPUT:  MSS - the subdomain matrix for the Schur system
%           MSE - the subdomain-edge matrix for the Schur system
%           MWW - the wirebasket matrix for the Schur system
%           MWE - the wirebasket-edge matrix for the Schur system
%           MEE - the edge matrix for the Schur system
%           LW - the wirebasket part of the Lagrangian multiplier
%           LE - the edge part of the Lagrangian multiplier
%           LW - the wirebasket part of the Lagrangian multiplier
%           LE - the edge part of the Lagrangian multiplier
%           b - the Schur vector
%           dsInd - the index for the subdomains
%           dwInd - the index for the wirebasket
%           deInd - the index for the edge
%           nV - mesh size
%           numDecompose - decomposition number
%           epsSchur - the convergence error of the conjugate method
%
%   OUTPUT: X - the solution
%           iterT - iteration time
%           t - running time


tic;

% %   compute solver for MSS and preconditioner
% MSSsolver = findsolverSS(MSS, numDecompose);
% presolver = findsolverPre(MWW, MWE, MEE_W);
%   compute XE
l = length(b);
XE = zeros(l, 1);
SX = SchurMultiplyWithSolver(MSEall, MWW, MWE, MEE, XE, MSSsolver, dssize, numDecompose);
% SX = SchurMultiplyWithSolver1(MSE, MWW, MWE, MEE, XE, MSSsolver, numDecompose);
rold = b - SX;
zold = Preconditioner(MWW, MWE, MEE_W, LW, LE, rold, presolver);
p = zold;
w = SchurMultiplyWithSolver(MSEall, MWW, MWE, MEE, p, MSSsolver, dssize, numDecompose);
% w = SchurMultiplyWithSolver1(MSE, MWW, MWE, MEE, p, MSSsolver, numDecompose);
alpha = (rold' * zold) / (p' * w);
XE = XE + alpha * p;
rnew = rold - alpha * w;
iterT = 1;

while norm(rnew, 2) > epsSchur
    znew = Preconditioner(MWW, MWE, MEE_W, LW, LE, rnew, presolver);
    beta = (rnew' * znew) / (rold' * zold);
    p = znew + beta * p;
    w = SchurMultiplyWithSolver(MSEall, MWW, MWE, MEE, p, MSSsolver, dssize, numDecompose);
%     w = SchurMultiplyWithSolver1(MSE, MWW, MWE, MEE, p, MSSsolver, numDecompose);
    alpha = (rnew' * znew) / (p' * w);
    XE = XE + alpha * p;
    rold = rnew;
    rnew = rnew - alpha * w;
    zold = znew;
    iterT = iterT + 1;
end

%   compute XS and XW
X0 = MSEall * XE;
for i = 1:numDecompose
    V0{i} = CS{i} - X0(dssize{i}, 1);
end
XS = MSSsolver.solve(V0);
XW = MWW \ (CW - MWE * XE);

%   assemble X
X = zeros(nV, 1);
for i = 1:numDecompose
    X(dsInd{i}) = XS{i};
end
X(dwInd) = XW;
X(deInd) = XE;


t = toc;


end