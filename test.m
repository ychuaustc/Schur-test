%%	construct sparse system MX = C with given decomposition
%%
%%	Input:
%%          M:              M matrix
%%          CT:             C vector (tuple)
%%          DS:             list of subdomain in the decomposition
%%          DB:             list of boundary in the decomposition
%%          DE:             edges in the decomposition
%%          DW:             wirebasket in the decomposition
%%          DWB:            boundary of wirebasket in the decomposition           
%%          numDecompose:   number of decompositions
%%	Output:
%%          MSS:            matrix block for subdomain * subdomain (tuple)
%%          MSE:            matrix block for subdomain * edge (tuple)
%%          MWW:            matrix block for wirebasket * wirebasket
%%          MWE:            matrix block for wirebasket * edge
%%          MEE:            matrix block for edge * edge
%%          CST:            matrix block for subdomain (2d tuple)
%%          CWT:            matrix block for wirebasket (tuple)
%%          S:              Schur matrix
%%          bT:             Schur vector (tuple)
%%          dsInd:          index for subdomain (tuple)
%%          dwInd:          index for wirebasket (tuple)
%%          deInd:          index for edge (tuple)

function [S] = test(M, DS, DB, DE, DW, DWB, numDecompose)

%% compute decomposition index
%  subdomain
for i = 1:numDecompose
    dsInd{i} = find(DS{i}); % subdomain index
end
%  wirebasket
dwInd = find(DW); % wirebasket index
%  edge
deInd = find(DE); % edge index
%  boundary
for i = 1:numDecompose
    dbInd{i} = find(DB{i}); % boundary index
end
dwbInd = find(DWB); % wirebasket boundary index

%%  define matrix blocks for M
for i = 1:numDecompose
    MSS{i} = M(dsInd{i}, dsInd{i}); % matrix block for subdomain * subdomain
end
for i = 1:numDecompose
    MSE{i} = M(dsInd{i}, deInd); % matrix block for subdomain * edge
end
MWW = M(dwInd, dwInd); % matrix block for wirebasket * wirebasket
MWE = M(dwInd, deInd); % matrix block for wirebasket * edge
MEE = M(deInd, deInd); % matrix block for edge * edge

%%  compute Schur matrix S and Schur vector b
%   Schur matrix S
S = MEE;
for i = 1:numDecompose
    S = S - MSE{i}' * (MSS{i} \ MSE{i});
end
S = S - MWE' * (MWW \ MWE);

%%
end