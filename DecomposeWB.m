%%  
%%
%%  Input:
%%          
%%  Output:
%%          

function [DV, DS, DB, DE, DW, DWB, numDS, numDecompose] = DecomposeWB(MC, B, DE, nV, nv)

%%
DW = DE;
DW(B) = 0;              %
dwInd = find(DW)';      %
dwbIndTemp = find(sum(MC(:, dwInd), 2))';
DWB = zeros(nV, 1);
DWB(dwbIndTemp) = 1;
DWB(dwInd) = 0;         %
dwbInd = find(DWB)';    %

%%
M = MC;
M(dwbInd, :) = 0;
M(:, dwbInd) = 0;
BB = zeros(nV, 1);
BB(B) = 1;
bwbInd = find(BB & DWB);
M(bwbInd, bwbInd) = 0;

%%
NewB = zeros(nV, 1);
NewB(B) = 1;
NewB(dwbInd) = 1;
qInd = find(sum(M, 2))';
NewB(qInd) = 0;
newbInd = find(NewB);

%%
[DV, DS, DB, DE, ~, ~, numDS, numDecompose] = DecomposeUT(nV, nv, newbInd, M);

%%
end