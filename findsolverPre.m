function [presolver] = findsolverPre(MWW, MWE, MEE_W)
M0{1} = [MWW MWE;MWE' MEE_W];
nonzero_ele{1} = nonzeros(M0{1});
presolver = batch_splsolver(M0, 'lu');
presolver.refactorize(nonzero_ele);