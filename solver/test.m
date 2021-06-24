% test script
A = cell(3, 1);
b = cell(3, 1);
x_ref = cell(3,1);
nonzero_ele = cell(3, 1);

A{1} = sparse([1 0 0 0 1
               0 1 0 1 0
               0 0 1 0 0
               1 0 0 1 0
               0 0 0 0 1]);
A{2} = sparse([1/2 1/2 0 0
               1/2 1/2 1/2 0
               0 1/2 1/2 1/2
               0 0 1/2 1/2]);
A{3} = sparse([1 2 3 0
               0 2 2 0
               0 0 1 0
               0 0 0 1]);

for i = 1:3
    b{i} = rand(size(A{i}, 1), 1);
    nonzero_ele{i} = nonzeros(A{i});
    x_ref{i} = A{i}\b{i};
end

% %% test full solve
% fprintf('test for full solve:\n');
% x = batch_splsolver.fullsolve(A, b, 'lu');
% err = 0;
% for i = 1:3
%     err = err + (norm(x{i} - x_ref{i}))^2;
% end
% fprintf('error = %f\n', sqrt(err));

%% create solver
solver = batch_splsolver(A, 'lu');
% %% test numerical factor & solve
% [x, f] = solver.refactor_solve(nonzero_ele, b);
% fprintf('test for factor & solve:\n');
% err = 0;
% for i = 1:3
%     err = err + (norm(x{i} - x_ref{i}))^2;
% end
% fprintf('error = %f\n', sqrt(err));
%% test refactorize
% for k = 1:3
%     nonzero_ele{k} = rand(length(nonzero_ele{k}), 1);
%     [i, j , v] = find(A{k});
%     A{k} = sparse(i, j, nonzero_ele{k});
%     x_ref{k} = A{k}\b{k};
% end
solver.refactorize(nonzero_ele);
%% resolve
fprintf('test for resolve:\n');
x = solver.solve(b);
err = 0;
for i = 1:3
    err = err + (norm(x{i} - x_ref{i}))^2;
end
fprintf('error = %f\n', sqrt(err));
