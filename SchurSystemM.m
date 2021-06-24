function [MSS, MSE, MSEall, MWW, MWE, MEE] = SchurSystemM(M, dsInd, dwInd, deInd, dsIndall, numDecompose)
%
%   this function computes the index and matrix blocks for the Schur system upon a given
%   decomposition
%
%   INPUT:  M - the coefficient matrix
%           dsInd - the index for the subdomains
%           dwInd - the index for the wirebasket
%           deInd - the index for the edge
%
%   OUTPUT: MSS - the subdomain matrix for the Schur system
%           MSE - the subdomain-edge matrix for the Schur system
%           MWW - the wirebasket matrix for the Schur system
%           MWE - the wirebasket-edge matrix for the Schur system
%           MEE - the edge matrix for the Schur system


%   compute the matrix blocks
for i = 1:numDecompose
    MSS{i} = M(dsInd{i}, dsInd{i});	% the matrix blocks for subdomains
end
for i = 1:numDecompose
    MSE{i} = M(dsInd{i}, deInd);	% matrix block for subdomain-edge
end
MSEall = M(dsIndall, deInd);
MWW = M(dwInd, dwInd);	% the matrix block for the wirebasket
MWE = M(dwInd, deInd);	% the matrix block for wirebasket-edge
MEE = M(deInd, deInd);	% the matrix block for the edge


end