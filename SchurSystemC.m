function [CS, CW, b] = SchurSystemC(MSS, MSE, MWW, MWE, C, dsInd, dwInd, deInd, numDecompose)
%
%   this function computes the right hand side of the Schur system
%
%   INPUT:  MSS - the subdomain matrix for the Schur system
%           MSE - the subdomain-edge matrix for the Schur system
%           MWW - the wirebasket matrix for the Schur system
%           MEE - the edge matrix for the Schur system
%           C - the right hand side of the system
%           dsInd - the index for the subdomains
%           dwInd - the index for the wirebasket
%           deInd - the index for the edge
%           numDecompose - decomposition number
%
%   OUTPUT: CS - the subdomain part of the right hand side
%           CW - the wirebasket part of the right hand side
%           b - the Schur vector


for i = 1:numDecompose
	CS{i} = C(dsInd{i});
end
CW = C(dwInd);
CE = C(deInd);

%   compute the Schur vector b
b = CE;
for i = 1:numDecompose
	b = b - MSE{i}' * (MSS{i} \ CS{i});
end
b = b - MWE' * (MWW \ CW);


end