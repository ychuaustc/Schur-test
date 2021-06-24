function [CS, CW, b] = SchurSystemC1(MSEall, MWW, MWE, C, dsInd, dwInd, deInd, solver, solverw, numDecompose)
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
X = solver.solve(CS);
Xall = cell2mat(X');
b = b - MSEall' * Xall;
% V0{1} = CW;
% U0 = solverw.solve(V0);
% b = b - MWE' * U0{1};
b = b - MWE' * (MWW \ CW);


end