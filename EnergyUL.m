function [EnUL] = EnergyUL(FCotTheta, FLVX, FLVY, FXU, FYU)
%
%   this function computes the ARAP energy with V and U
%
%   INPUT:  FCotTheta - the cotangent value of the face matrix
%           FLVX - assign the x-coordinate value of L(V) to the face matrix
%           FLVY - assign the y-coordinate value of L(V) to the face matrix
%           FXU - x-coordinates for the face matrix-U
%           FYU - y-coordinates for the face matrix-U
%
%   OUTPUT: EnUL - the ARAP energy E(U, L)


EnUL = 0;

for i = 1:3
    ip1 = mod(i, 3) + 1;	% (i + 1) mod 3

    CXi = (FXU(:, i) - FXU(:, ip1)) - (FLVX(:, i) - FLVX(:, ip1));
    CYi = (FYU(:, i) - FYU(:, ip1)) - (FLVY(:, i) - FLVY(:, ip1));
    Ei = 0.5 * sum(FCotTheta(:, i) .* (CXi .* CXi + CYi .* CYi));
    EnUL = EnUL + Ei;
end


end