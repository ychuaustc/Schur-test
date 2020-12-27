%%  return a decomposition of given mesh
%%
%%  Input:
%%  Output:

function [DV, DS, DE, DW, DWB, numDecompose] = DecomposeTest2(nV)

%%
sqnV = sqrt(nV);
%
v1 = [1:sqnV:(sqnV - 1) * sqnV + 1];
DW = zeros(nV, 1);
DW(v1) = 1;
%
v2 = v1 + 1;
DWB = zeros(nV, 1);
DWB(v2) = 1;
DE = zeros(nV, 1);
DE(v2) = 1;
%
DS{1} = ones(nV, 1);
DS{1}(v1) = 0;
DS{1}(v2) = 0;
%
DV{1} = DS{1};
DV{1}(v2) = 1;
%
numDecompose = 1;

%%
end