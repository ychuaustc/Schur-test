%%  
%%
%%  Input:
%%          
%%  Output:
%%          

function [DS, DB, DE, DW, DWB] = ModifyDecompose(DV, numDecompose, MC, nV)

%%
for i = 1:numDecompose
    DS{i} = ones(nV, 1);
end
DW = zeros(nV, 1);
DWB = zeros(nV, 1);
DE = zeros(nV, 1);
for i = 1:numDecompose
    for j = 1:numDecompose
        if i ~= j
            dvInd = find(DV{j})';
            DS{i}(dvInd) = 0;
        end
    end
    DB{i} = double(DV{i} > DS{i});
end
for i = 1:numDecompose
    dbInd = find(DB{i})';
    DW(dbInd) = 1;                      %
end
dwInd = find(DW)';
DWB(find(sum(MC(:, dwInd), 2))') = 1;
DWB(dwInd) = 0;                         %
%
sqnV = sqrt(nV);
id1 = ((sqnV + 1) / 2) * sqnV + (sqnV - 1) / 2;
id2 = ((sqnV - 1) / 2) * sqnV - (sqnV - 1) / 2 + 1;
DW(id1) = 1;
DW(id2) = 1;
DWB(id1) = 0;
DWB(id2) = 0;
%
DE = DWB;                               %
DWWB = double(DW | DWB);
for i = 1:numDecompose
    DS{i} = double(DS{i} > DWWB);
end                                     %

%%
end