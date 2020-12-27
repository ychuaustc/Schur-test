%%  return a decomposition of given mesh
%%
%%  Input:
%%  Output:

function [DV, DS, DE, DW, DWB, numDecompose] = DecomposeTest(nV)

%%
sqnV = sqrt(nV);
%
dwInd = 3:sqnV:nV - sqnV + 3;
dwbInd1 = 2:sqnV:nV - sqnV + 2;
dwbInd2 = 4:sqnV:nV - sqnV + 4;
dsInd1 = 1:sqnV:nV - sqnV + 1;
dsInd2 = setdiff(1:nV, dwInd);
dsInd2 = setdiff(dsInd2, dwbInd1);
dsInd2 = setdiff(dsInd2, dwbInd2);
dsInd2 = setdiff(dsInd2, dsInd1);
%
DW = zeros(nV, 1);
DV{1} = zeros(nV, 1);
DV{2} = zeros(nV, 1);
DS{1} = zeros(nV, 1);
DS{2} = zeros(nV, 1);
DWB = zeros(nV, 1);
DE = zeros(nV, 1);
%
DW(dwInd) = 1;
DWB(dwbInd1) = 1;
DWB(dwbInd2) = 1;
DE = DWB;
DS{1}(dsInd1) = 1;
DS{2}(dsInd2) = 1;
DV{1} = DS{1};
DV{1}(dwbInd1) = 1;
DV{2} = DS{2};
DV{2}(dwbInd2) = 1;
numDecompose = 2;






% %%
% sqnV = sqrt(nV);
% %
% dwInd1 = 1:sqnV:nV - sqnV + 1;
% dwInd2 = 5:sqnV:nV - sqnV + 5;
% dwInd3 = sqnV - 4:sqnV:nV - sqnV + sqnV - 4;
% dwInd4 = sqnV:sqnV:nV - sqnV + sqnV;
% dwbInd1 = 2:sqnV:nV - sqnV + 2;
% dwbInd2 = 4:sqnV:nV - sqnV + 4;
% dwbInd3 = 6:sqnV:nV - sqnV + 6;
% dwbInd4 = sqnV - 5:sqnV:nV - sqnV + sqnV - 5;
% dwbInd5 = sqnV - 3:sqnV:nV - sqnV + sqnV - 3;
% dwbInd6 = sqnV - 1:sqnV:nV - sqnV + sqnV - 1;
% %
% DW = zeros(nV, 1);
% DWB = zeros(nV, 1);
% DE = zeros(nV, 1);
% DV{1} = zeros(nV, 1);
% DV{2} = zeros(nV, 1);
% DV{3} = zeros(nV, 1);
% DS{1} = zeros(nV, 1);
% DS{2} = zeros(nV, 1);
% DS{3} = zeros(nV, 1);
% %
% DW(dwInd1) = 1;
% DW(dwInd2) = 1;
% DW(dwInd3) = 1;
% DW(dwInd4) = 1;
% DWB(dwbInd1) = 1;
% DWB(dwbInd2) = 1;
% DWB(dwbInd3) = 1;
% DWB(dwbInd4) = 1;
% DWB(dwbInd5) = 1;
% DWB(dwbInd6) = 1;
% DE = DWB;
% DS{1}(3:sqnV:nV - sqnV + 3) = 1;
% for i = 7:sqnV - 6
%     DS{2}(i:sqnV:nV - sqnV + i) = 1;
% end
% DS{3}(sqnV - 2:sqnV:nV - sqnV + sqnV - 2) = 1;
% DV{1} = DS{1};
% DV{2} = DS{2};
% DV{3} = DS{3};
% DV{1}(dwbInd1) = 1;
% DV{1}(dwbInd2) = 1;
% DV{2}(dwbInd3) = 1;
% DV{2}(dwbInd4) = 1;
% DV{3}(dwbInd5) = 1;
% DV{3}(dwbInd6) = 1;
% numDecompose = 3;

%%
end