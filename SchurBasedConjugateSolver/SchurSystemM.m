function [MSS, MSE, MWW, MWE, MEE, dsInd, dwInd, deInd] = SchurSystemM(M, DS, DE, DW, numDecompose)
%
%   this function computes the index and matrix blocks for the Schur system upon a given
%   decomposition
%
%   INPUT:  M - the coefficient matrix
%           DS - subdomain list
%           DE - edge of all decompositions
%           DW - wirebasket set
%           numDecompose - decomposition number
%
%   OUTPUT: MSS - the subdomain matrix for the Schur system
%           MSE - the subdomain-edge matrix for the Schur system
%           MWW - the wirebasket matrix for the Schur system
%           MWE - the wirebasket-edge matrix for the Schur system
%           MEE - the edge matrix for the Schur system
%           dsInd - the index for the subdomains
%           dwInd - the index for the wirebasket
%           deInd - the index for the edge


%   the index in the decomposition
for i = 1:numDecompose
    dsInd{i} = find(DS{i});	% the index for the subdomains
end
dwInd = find(DW);	% the index for the wirebasket
deInd = find(DE);	% the index for the edge


%   compute the matrix blocks
for i = 1:numDecompose
    MSS{i} = M(dsInd{i}, dsInd{i});	% the matrix blocks for subdomains
end
for i = 1:numDecompose
    MSE{i} = M(dsInd{i}, deInd);	% matrix block for subdomain-edge
end
MWW = M(dwInd, dwInd);	% the matrix block for the wirebasket
MWE = M(dwInd, deInd);	% the matrix block for wirebasket-edge
MEE = M(deInd, deInd);	% the matrix block for the edge


end