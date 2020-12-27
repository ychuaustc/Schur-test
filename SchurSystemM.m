%%	construct sparse system MX = C with given decomposition
%%
%%	Input:
%%          M:              M matrix
%%          DS:             list of subdomain in the decomposition
%%          DB:             list of boundary in the decomposition
%%          DE:             edges in the decomposition
%%          DW:             wirebasket in the decomposition
%%          DWB:            boundary of wirebasket in the decomposition           
%%          numDecompose:   number of decompositions
%%	Output:
%%          MSS:            matrix block for subdomain * subdomain
%%          MSE:            matrix block for subdomain * edge (tuple)
%%          MWW:            matrix block for wirebasket * wirebasket
%%          MWE:            matrix block for wirebasket * edge
%%          MEE:            matrix block for edge * edge
%%          dsInd:          index for subdomain
%%          dwInd:          index for wirebasket
%%          deInd:          index for edge

function [MSS, MSE, MWW, MWE, MEE, dsInd, dwInd, deInd] = SchurSystemM(M, DS, DE, DW, numDecompose)

%% compute decomposition index
%  subdomain
for i = 1:numDecompose
    dsInd{i} = find(DS{i}); % subdomain index
end
%  wirebasket
dwInd = find(DW); % wirebasket index
%  edge
deInd = find(DE); % edge index

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

%%
end