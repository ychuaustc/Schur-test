function [] = ShowResult1(iterARAP_D, t_D, iterARAP_CP, iterSchur_CP, t_CP)
%
%   this function show the result of the parameterization
%
%   INPUT:  iterARAP_D - local-global iteration time (direct solver)
%           t_D - running time per local-global iteration (direct solver)
%           iterARAP_C - local-global iteration time (conjugate gradient solver)
%           iterSchur_C - iteration time for Schur method per local-global iteration (conjugate gradient solver)
%           t_C - running time per local-global iteration (preconditioned conjugate gradient solver)
%           iterARAP_CP - local-global iteration time (preconditioned conjugate gradient solver)
%           iterSchur_CP - iteration time for Schur method per local-global iteration (preconditioned conjugate gradient 
%                          solver)
%           t_CP - running time per local-global iteration (preconditioned conjugate gradient solver)
%
%   OUTPUT: 


fprintf('The local-global iteration time is: direct solver: %d, PCG solver: %d \n', iterARAP_D, iterARAP_CP);
fprintf('The Schur iteration/running time per local-global step of all solvers are given as follows:\n');

C1{1, 1} = 'Local-Global Iterations';
C1{2, 1} = 'Iterations/Local-Global Step';
C1{3, 1} = 'Running Time/Local Global Step';
C2 = [iterARAP_D; 0; t_D];
C3 = [iterARAP_CP; iterSchur_CP; t_CP];
T = table(C2, C3, 'VariableNames', {'Direct Solver', 'PCG'});
T.Properties.RowNames = C1


end