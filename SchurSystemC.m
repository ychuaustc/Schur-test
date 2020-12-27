%%	construct sparse system MX = C with given decomposition
%%
%%	Input:
%%          C:              C vector
%%          dsInd:          index for subdomain
%%          dwInd:          index for wirebasket
%%          deInd:          index for edge
%%          numDecompose:   number of decompositions
%%	Output:
%%          CS:             matrix block for subdomain
%%          CW:             matrix block for wirebasket
%%          b:              Schur vector

function [CS, CW, b] = SchurSystemC(MSS, MSE, MWW, MWE, C, dsInd, dwInd, deInd, numDecompose)

%%  define matrix blocks for CT
for i = 1:numDecompose
	CS{i} = C(dsInd{i});
end
CW = C(dwInd);
CE = C(deInd);

%%	Schur vector b
b = CE;
for i = 1:numDecompose
	b = b - MSE{i}' * (MSS{i} \ CS{i});
end
b = b - MWE' * (MWW \ CW);

%%
end