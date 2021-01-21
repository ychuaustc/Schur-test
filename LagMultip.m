function [LW, LE] = LagMultip(FX, FY, DW, DWB, dwInd, deInd, Face, nV, nF)
%
%   this function computes the lagrangian multiplier if the Schur system is singular
%
%   INPUT:  FX - x-coordinates for the face matrix
%           FY - y-coordinates for the face matrix
%           DW - wirebasket set
%           DWB - wirebasket boundary
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
IntegralPhi = 2 * (1 / 3) * 1 * AreaF;  % list of the integral of piecewise linear functions (P1-element) on each face

FInWWB1 = zeros(nF, 1);
FInWWB2 = zeros(nF, 1);
FInWWB3 = zeros(nF, 1);
Face1 = Face(:, 1);
Face2 = Face(:, 2);
Face3 = Face(:, 3);
dwwbInd = find(DW | DWB)';
for i = dwwbInd
    finwwb1 = find(Face1 == i)';
    FInWWB1(finwwb1) = 1;
    finwwb2 = find(Face2 == i)';
    FInWWB2(finwwb2) = 1;
    finwwb3 = find(Face3 == i)';
    FInWWB3(finwwb3) = 1;
end
finwInd = find(FInWWB1 & FInWWB2 & FInWWB3)';   % faces which has points in the wirebasket set

L = zeros(nV, 1);
for i = finwInd
    for j = Face(i, :)
        L(j) = L(j) + IntegralPhi(i);
    end
end
LW = L(dwInd);
LE = L(deInd);


end