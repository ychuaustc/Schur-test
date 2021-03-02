function [DS, DE, DW, numDecompose] = Decomp1(nV)

numDecompose = 1;

nv = sqrt(nV);

DS{1} = ones(nV, 1);
DW = zeros(nV, 1);
DE = zeros(nV, 1);

DW(1:nv:(nv - 1) * nv + 1) = 1;
DE(2:nv:(nv - 1) * nv  +2) = 1;
dw = find(DW)';
de = find(DE)';
DS{1}(dw) = 0;
DS{1}(de) = 0;

end