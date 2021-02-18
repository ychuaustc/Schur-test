function [DV, DS, DB, DE, DW, DWB, numDecompose] = DecomposeTest(nV)


numDecompose = 1;
DS{1} = ones(nV, 1);
DB{1} = zeros(nV, 1);
DV{1} = zeros(nV, 1);
DE = zeros(nV, 1);
DW = zeros(nV, 1);
DWB = zeros(nV, 1);

nv = sqrt(nV);
w1 = 1:nv:(nv - 1) * nv + 1;
% w1 = [];
w2 = nv:nv:nV;
% w2 = [];
wb1 = 2:nv:(nv - 1) * nv + 2;
% wb1 = [];
wb2 = nv - 1:nv:nV - 1;
% wb2 = [];

DW(w1) = 1;
DW(w2) = 1;
DWB(wb1) = 1;
DWB(wb2) = 1;
DE = DWB;
DS{1}(w1) = 0;
DS{1}(w2) = 0;
DS{1}(wb1) = 0;
DS{1}(wb2) = 0;
DV{1} = DS{1};
DV{1}(wb1) = 1;
DV{1}(wb2) = 1;


end