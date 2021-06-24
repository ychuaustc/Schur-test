function [LW, LE] = LagMultip(FX, FY, dwInd, deInd, finweInd, Face, nV)
%
%   this function computes the lagrangian multiplier if the Schur system is singular
%
%   INPUT:  FX - x-coordinates for the face matrix
%           FY - y-coordinates for the face matrix
%           DW - wirebasket set
%           DE - edge of all decomposition
%           dwInd - the index for the wirebasket
%           deInd - the index for the edge
%           Face - mesh faces
%           nV - mesh size
%           nF - number of faces
%
%   OUTPUT: LW - the wirebasket part of the Lagrangian multiplier
%           LE - the edge part of the Lagrangian multiplier


F1X = FX(:, 1);
F2X = FX(:, 2);
F3X = FX(:, 3);
F1Y = FY(:, 1);
F2Y = FY(:, 2);
F3Y = FY(:, 3);

E2X = F2X - F1X;
E2Y = F2Y - F1Y;
E3X = F3X - F1X;
E3Y = F3Y - F1Y;
AreaF = 0.5 * abs(E2X .* E3Y - E3X .* E2Y);
IntegralPhi = 2 * (1 / 3) * 1 * AreaF;  % list of the integral of P1-element on each face

L = zeros(nV, 1);
for i = finweInd
    for j = Face(i, :)
        L(j) = L(j) + IntegralPhi(i);
    end
end
LW = L(dwInd);
LE = L(deInd);


end