function [X, iterT, t] = SchurConjSolver(MSS, MSE, MWW, MWE, MEE, CS, CW, b, solverMSS, dsInd, dwInd, deInd, nV, numDecompose, epsSchur)
%
%   this function solves the Schur system using conjugate gradient method
%
%   INPUT:  MSS - the subdomain matrix for the Schur system
%           MSE - the subdomain-edge matrix for the Schur system
%           MWW - the wirebasket matrix for the Schur system
%           MWE - the wirebasket-edge matrix for the Schur system
%           MEE - the edge matrix for the Schur system
%           CS - the subdomain part of the right hand side
%           CW - the wirebasket part of the right hand side
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
Sx = SchurMultiply(MSS, MSE, MWW, MWE, MEE, XE, solverMSS, numDecompose);
r = b - Sx;
p = r;
rsold = r' * r;

iterT = 0;
for i = 1:l
    Sp = SchurMultiply(MSS, MSE, MWW, MWE, MEE, p, solverMSS, numDecompose);
    alpha = rsold / (p' * Sp);
    XE = XE + alpha * p;
    r = r - alpha * Sp;
    rsnew = r' * r;
    iterT = iterT + 1;
    if sqrt(rsnew) < epsSchur
        break;
    end
    p = r + (rsnew / rsold) * p;
    rsold = rsnew;
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