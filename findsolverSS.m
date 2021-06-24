function [MSSsolver] = findsolverSS(MSS, numDecompose)
for i = 1:numDecompose
    nonzero_ele{i} = nonzeros(MSS{i});
end
MSSsolver = batch_splsolver(MSS, 'lu');
MSSsolver.refactorize(nonzero_ele);