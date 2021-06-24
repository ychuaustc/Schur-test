function [prehsolver] = findsolverPreH(M)
M0{1} = M;
nonzero_ele{1} = nonzeros(M0{1});
prehsolver = batch_splsolver(M0, 'lu');
prehsolver.refactorize(nonzero_ele);