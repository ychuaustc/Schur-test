function [presolver] = findsolverW(MWW)
M0{1} = MWW;
nonzero_ele{1} = nonzeros(M0{1});
presolver = batch_splsolver(M0, 'lu');
presolver.refactorize(nonzero_ele);