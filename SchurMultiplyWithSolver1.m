function U = SchurMultiplyWithSolver1(MSE, MWW, MWE, MEE, V, solver, numDecompose)
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


% U = MEE * V;
% for i = 1:numDecompose
%     U = U - MSE{i}' * (MSS{i} \ (MSE{i} * V));
% %     U = U - MSE{i}' * solverMSS{i}.refactor_solve(MSS{i}, MSE{i} * V);
% end
% U = U - MWE' * (MWW \ (MWE * V));


U = MEE * V;
for i = 1:numDecompose
    b{i} = MSE{i} * V;
end
X = solver.solve(b);
for i = 1:numDecompose
    U = U - MSE{i}' * X{i};
end
U = U - MWE' * (MWW \ (MWE * V));


end