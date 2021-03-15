function U = SchurMultiply(MSS, MSE, MWW, MWE, MEE, V, solverMSS, numDecompose)
%
%   this function computes the Schur matrix S times the vector V
%
%   INPUT:  MSS - the subdomain matrix for the Schur system
%           MSE - the subdomain-edge matrix for the Schur system
%           MWW - the wirebasket matrix for the Schur system
%           MWE - the wirebasket-edge matrix for the Schur system
%           MEE - the edge matrix for the Schur system
%           V - the vector
%           numDecompose - decomposition number
%
%   OUTPUT: U - S times V


U = MEE * V;
for i = 1:numDecompose
% parfor i = 1:numDecompose
%     U = U - MSE{i}' * (MSS{i} \ (MSE{i} * V));
    U = U - MSE{i}' * solverMSS{i}.refactor_solve(MSS{i}, MSE{i} * V);
end
U = U - MWE' * (MWW \ (MWE * V));


end