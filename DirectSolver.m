%%  solve the Schur preconditioner problem using a direct solver (single core)
%%
%%  Input:  MS: subdomain matrix block of the coefficient matrix (tuple of sparse matrix)
%%          MSE: subdomain intersect edges matrix block of the coefficient matrix (tuple of sparse matrix)
%%          MW: wirebasket matrix block of the coefficient matrix (sparse matrix)
%%          MWE: wirebasket intersect edges matrix block of the coefficient matrix (sparse matrix)
%%          ME: edges matrix block of the coefficient matrix (sparse matrix)
%%          CS: subdomain vector for the RHS of the linear system (tuple of column vector)
%%          CW: wirebasket vector for the RHS of the linear system (column vector)
%%          CE: edges vector for the RHS of the linear system (column vector)
%%          numDecompose: number of decomposition (constant)
%%  Output: XS: subdomain part of the solution(tuple of colume vector)
%%          XW: wirebasket part of the solution (column vector)
%%          XE: edges part of the solution (column vector)
%%          t: running time

function [XS, XW, XE, t] = DirectSolver(MS, MSE, MW, MWE, ME, CS, CW, CE, numDecompose)

%%
tic;

%% compute XE
%
b = CE;
for i = 1:numDecompose
    b = b - MSE{i}' * (MS{i} \ CS{i});
end
b = b - MWE' * (MW \ CW);
%
S = ME;
for i = 1:numDecompose
    S = S - MSE{i}' * (MS{i} \ MSE{i});
end
S = S - MWE' * (MW \ MWE);
%
XE = S \ b;

%% compute XS with XI
for i = 1:numDecompose
    XS{i} = MS{i} \ (CS{i} - MSE{i} * XE);
end
if sum(MW) ~= 0
    XW = MW \ (CW - MWE * XE);
else
    XW = 0;
end

%%
t = toc;

end