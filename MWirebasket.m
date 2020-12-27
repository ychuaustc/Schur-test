%%  
%%
%%  Input:
%%          
%%  Output:
%%          

function [MEE_W] = MWirebasket(FCotTheta, Face, DW, DWB, deInd, nV, nF)

%%
%
FNotInWWB = zeros(nF, 1);
Face1 = Face(:, 1);
Face2 = Face(:, 2);
Face3 = Face(:, 3);
dwwbInd = find(DW | DWB)';
notindwwbInd = setdiff(1:nV, dwwbInd);
for i = notindwwbInd
    fnotinwwb1 = find(Face1 == i)';
    FNotInWWB(fnotinwwb1) = 1;
    fnotinwwb2 = find(Face2 == i)';
    FNotInWWB(fnotinwwb2) = 1;
    fnotinwwb3 = find(Face3 == i)';
    FNotInWWB(fnotinwwb3) = 1;
end
fnotinwInd = find(FNotInWWB)';
%
FT = FCotTheta;
FT(fnotinwInd, :) = 0;
MT = sparse(Face, Face(:, [2 3 1]), FT(:, [3 1 2]), nV, nV);
MTemp = sparse(1:nV, 1:nV, sum(MT + MT', 2), nV, nV) - (MT + MT');
MEE_W = MTemp(deInd, deInd);

%%
end