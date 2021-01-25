function [X, iterT, t] = SchurConjPreSolver(MSS, MSE, MWW, MWE, MEE, MEE_W, CS, CW, b, dsInd, dwInd, deInd, nV, numDecompose, epsSchur)
%
%   this function solves the Schur system using a preconditioned conjugate method
%
%   INPUT:  MSS - the subdomain matrix for the Schur system
%           MSE - the subdomain-edge matrix for the Schur system
%           MWW - the wirebasket matrix for the Schur system
%           MWE - the wirebasket-edge matrix for the Schur system
%           MEE - the edge matrix for the Schur system
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

%   compute XE
l = length(b);
XE = zeros(l, 1);
SX = SchurMultiply(MSS, MSE, MWW, MWE, MEE, XE, numDecompose);
rold = b - SX;
zold = Preconditioner(MWW, MWE, MEE_W, rold);
p = zold;
w = SchurMultiply(MSS, MSE, MWW, MWE, MEE, p, numDecompose);
alpha = (rold' * zold) / (p' * w);
XE = XE + alpha * p;
rnew = rold - alpha * w;
iterT = 1;

while norm(rnew, 2) > epsSchur
    znew = Preconditioner(MWW, MWE, MEE_W, rnew);
    beta = (rnew' * znew) / (rold' * zold);
    p = znew + beta * p;
    w = SchurMultiply(MSS, MSE, MWW, MWE, MEE, p, numDecompose);
    alpha = (rnew' * znew) / (p' * w);
    XE = XE + alpha * p;
    rold = rnew;
    rnew = rnew - alpha * w;
    zold = znew;
    iterT = iterT + 1;
end

%   compute XS and XW
for i = 1:numDecompose
    XS{i} = MSS{i} \ (CS{i} - MSE{i} * XE);
end
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