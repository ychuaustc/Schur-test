%%	solving sparse system MX = C via Schur compliment method (SX' = b) on decomposed mesh (conjugate version)
%%
%%	Input:
%%          MSS:            matrix block for subdomain * subdomain
%%          MSE:            matrix block for subdomain * edge
%%          MWW:            matrix block for wirebasket * wirebasket
%%          MWE:            matrix block for wirebasket * edge
%%          MEE:            matrix block for edge * edge
%%          CS:             matrix block for subdomain (tuple)
%%          CW:             matrix block for wirebasket
%%          b:              Schur vector
%%          numDecompose:   number of decomposition
%%          epsSchur:   value for Schur solver convergence check
%%
%%	Output:
%%          X:              X
%%          iterT:          iteration time
%%          t:              running time
%%          

function [X, iterT, t] = SchurConjSolver(MSS, MSE, MWW, MWE, MEE, CS, CW, b, dsInd, dwInd, deInd, nV, numDecompose, epsSchur)

%%
tic;

%%	compute XE
l = length(b);
XE = zeros(l, 1);
Sx = SchurMultiply(MSS, MSE, MWW, MWE, MEE, XE, numDecompose);
r = b - Sx;
p = r;
rsold = r' * r;
%
iterT = 0;
for i = 1:l
    Sp = SchurMultiply(MSS, MSE, MWW, MWE, MEE, p, numDecompose);
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

%%	compute XS and XW
for i = 1:numDecompose
    XS{i} = MSS{i} \ (CS{i} - MSE{i} * XE);
end
XW = MWW \ (CW - MWE * XE);

%%  assemble X
X = zeros(nV, 1);
for i = 1:numDecompose
    X(dsInd{i}) = XS{i};
end
X(dwInd) = XW;
X(deInd) = XE;

%%
t = toc;

%%
end