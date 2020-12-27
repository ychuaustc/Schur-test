%%  return a decomposition of given mesh
%%
%%  Input:
%%  Output:

function [DV, DS, DE, DW, DWB, numDecompose] = DecomposeTest1(MC, nV, sqnV)

%%  3   1
numDecompose = 2;
% 
%%
for i = 1:numDecompose
    DS{i} = zeros(nV, 1);
    DB{i} = zeros(nV, 1);
    DV{i} = zeros(nV, 1);
end
DV{numDecompose + 1} = zeros(nV, 1);
DE = zeros(nV, 1);
DW = zeros(nV, 1);
DWB = zeros(nV, 1);

%%
%
for i = 1:sqnV
    DS{1}([(i - 1) * sqnV + 4:(i - 1) * sqnV + (sqnV - 1) / 2]) = 1;
    DS{2}([(i - 1) * sqnV + (sqnV + 1) / 2 + 1:(i - 1) * sqnV + sqnV - 3]) = 1;
end
for i = 1:numDecompose
    dsInd{i} = find(DS{i});
    DV{i}(dsInd{i}) = 1;
end
%
for i = 1:numDecompose
    vecB = setdiff(find(sum(MC(:, dsInd{i}), 2))', dsInd{i});
    DB{i}(vecB) = 1;
    DE(vecB) = 1;
    DV{i}(vecB) = 1;
end
%
vecW1 = union([1:sqnV:(sqnV - 1) * sqnV + 1], [2:sqnV:(sqnV - 1) * sqnV + 2]);
vecW2 = union([sqnV - 1:sqnV:(sqnV - 1) * sqnV + sqnV - 1], [sqnV:sqnV:(sqnV - 1) * sqnV + sqnV]);
vecW = union(vecW1, vecW2);
DW(vecW) = 1;
vecWB = setdiff(find(sum(MC(:, vecW), 2))', vecW);
DWB(vecWB) = 1;
DE(vecWB) = 1;
DV{numDecompose + 1}(vecW) = 1;
DV{numDecompose + 1}(vecWB) = 1;

%%
end