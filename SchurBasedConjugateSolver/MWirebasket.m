function [MEE_W] = MWirebasket(FCotTheta, Face, DW, DE, deInd, nV, nF)
%
%   this function computes the wirebasket matrix block
%
%   INPUT:  FCotTheta - the cotangent value of the face matrix
%           Face - mesh faces
%           DW - wirebasket set
%           DE - edge of all decomposition
%           deInd - the index for the edge
%           nV - mesh size
%           nF - number of faces
%
%   OUTPUT: MEE_W - the wirebasket matrix block


FNotInWE = zeros(nF, 1);
Face1 = Face(:, 1);
Face2 = Face(:, 2);
Face3 = Face(:, 3);
dweInd = find(DW | DE)';
notinweInd = setdiff(1:nV, dweInd);
for i = notinweInd
    fnotinwe1 = find(Face1 == i)';
    FNotInWE(fnotinwe1) = 1;
    fnotinwe2 = find(Face2 == i)';
    FNotInWE(fnotinwe2) = 1;
    fnotinwe3 = find(Face3 == i)';
    FNotInWE(fnotinwe3) = 1;
end
fnotinweInd = find(FNotInWE)';

FT = FCotTheta;
FT(fnotinweInd, :) = 0;
MT = sparse(Face, Face(:, [2 3 1]), FT(:, [3 1 2]), nV, nV);
MTemp = sparse(1:nV, 1:nV, sum(MT + MT', 2), nV, nV) - (MT + MT');
MEE_W = MTemp(deInd, deInd);


end