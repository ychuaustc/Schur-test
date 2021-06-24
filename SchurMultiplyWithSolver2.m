function U = SchurMultiplyWithSolver2(MSEall, MWW, MWE, MEE, V, solver, solverw, dssize, numDecompose)
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
b0 = MSEall * V;
for i = 1:numDecompose
    b{i} = b0(dssize{i}, 1);
end
X = solver.solve(b);
% Xall = [];
% for i = 1:numDecompose
%     Xall = [Xall; X{i}];
% end
Xall = cell2mat(X');
U = U - MSEall' * Xall;
% V0{1} = MWE * V;
% U0 = solverw.solve(V0);
% U = U - MWE' * U0{1};
U = U - MWE' * (MWW \ (MWE * V));


end