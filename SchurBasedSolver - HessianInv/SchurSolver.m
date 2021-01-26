%%	solving sparse system MX = C via Schur compliment method (SX' = b) on decomposed mesh
%%
%%	Input:
%%          MSS:            matrix block for subdomain * subdomain
%%          MSE:            matrix block for subdomain * edge
%%          MWW:            matrix block for wirebasket * wirebasket
%%          MWE:            matrix block for wirebasket * edge
%%          CS:             matrix block for subdomain (tuple)
%%          CW:             matrix block for wirebasket
%%          CE:             matrix block for edge
%%          S:              Schur matrix
%%          b:              Schur vector
%%          dsInd:          index for subdomain (tuple)
%%          dwInd:          index for wirebasket (tuple)
%%          deInd:          index for edge (tuple)
%%          nV:             mesh size
%%          numDecompose:   number of decompositions
%%
%%	Output:
%%          X:      X
%%          t:      running time
%%          

function [X, t] = SchurSolver(MSS, MSE, MWW, MWE, CS, CW, S, b, dsInd, dwInd, deInd, nV, numDecompose)

%%
tic;

%%	compute XE
XE = S \ b;

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