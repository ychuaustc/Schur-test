function [MEE_W] = MWirebasket(FCotTheta, Face, DW, DWB, deInd, nV, nF)
%
%   this function computes the wirebasket matrix block
%
%   INPUT:  FCotTheta - the cotangent value of the face matrix
%           Face - mesh faces
%           DW - wirebasket set
%           DWB - wirebasket boundary
%           deInd - the index for the edge
%           nV - mesh size
%           nF - number of faces
%
%   OUTPUT: MEE_W - the wirebasket matrix block


FNotInWWB = zeros(nF, 1);
Face1 = Face(:, 1);
Face2 = Face(:, 2);
Face3 = Face(:, 3);
dwwbInd = find(DW | DWB)';
notinwwbInd = setdiff(1:nV, dwwbInd);
for i = notinwwbInd
    fnotinwwb1 = find(Face1 == i)';
    FNotInWWB(fnotinwwb1) = 1;
    fnotinwwb2 = find(Face2 == i)';
    FNotInWWB(fnotinwwb2) = 1;
    fnotinwwb3 = find(Face3 == i)';
    FNotInWWB(fnotinwwb3) = 1;
end
fnotinwInd = find(FNotInWWB)';

FT = FCotTheta;
FT(fnotinwInd, :) = 0;
MT = sparse(Face, Face(:, [2 3 1]), FT(:, [3 1 2]), nV, nV);
MTemp = sparse(1:nV, 1:nV, sum(MT + MT', 2), nV, nV) - (MT + MT');
MEE_W = MTemp(deInd, deInd);


end