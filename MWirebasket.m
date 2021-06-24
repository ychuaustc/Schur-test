function [MEE_W] = MWirebasket(FCotTheta, Face, deInd, finweInd, nV, nF)
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


fnotinweInd = setdiff(1:nF, finweInd)';

FT = FCotTheta;
FT(fnotinweInd, :) = 0;
MT = sparse(Face, Face(:, [2 3 1]), FT(:, [3 1 2]), nV, nV);
MTemp = sparse(1:nV, 1:nV, sum(MT + MT', 2), nV, nV) - (MT + MT');
MEE_W = MTemp(deInd, deInd);


end