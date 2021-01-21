function [FLVX, FLVY, MThetaLVX, MThetaLVY] = ArapL(FXV, FYV, FXU, FYU, FCotTheta, Face, nV, nF)
%
%   this function computes the rotation matrix for the local phase (V -> U)
%
%   INPUT:  FXV - x-coordinates for the face matrix-V
%           FYV - y-coordinates for the face matrix-V
%           FXU - x-coordinates for the face matrix-U
%           FYU - y-coordinates for the face matrix-U
%           FCotTheta - the cotangent value of the face matrix
%           FCotTheta - the cotangent value of the face matrix
%           nV - mesh size
%           nF - number of faces
%
%   OUTPUT: FLVX - assign the x-coordinate value of L(V) to the face matrix
%           FLVY - assign the y-coordinate value of L(V) to the face matrix
%           MThetaLVX - assign the x-coordinate value of cot(theta) * L(V) to the adjacent matrix
%           MThetaLVY - assign the y-coordinate value of cot(theta) * L(V) to the adjacent matrix


%   assemble the diagonal matrix with cross covariance matrix and do SVD
%   (a b)^T(c d) = (ac ad;bc bd) := (A B;C D)
A = zeros(nF, 1);
B = zeros(nF, 1);
C = zeros(nF, 1);
D = zeros(nF, 1);
for i = 1:3
    ip1 = mod(i, 3) + 1;	% (i + 1) mod 3

    FXU_imip1_i = FXU(:, i) - FXU(:, ip1);
    FYU_imip1_i = FYU(:, i) - FYU(:, ip1);
    ThetaFXU_imip1_i = FCotTheta(:, i) .* FXU_imip1_i;	% a
    ThetaFYU_imip1_i = FCotTheta(:, i) .* FYU_imip1_i;  % b
    FXV_imip1_i = FXV(:, i) - FXV(:, ip1);              % c
    FYV_imip1_i = FYV(:, i) - FYV(:, ip1);              % d
    
    % the matrix is of the form (a b;c d)
    A = A + ThetaFXU_imip1_i .* FXV_imip1_i;	% A = ac
    B = B + ThetaFXU_imip1_i .* FYV_imip1_i;	% B = ad
    C = C + ThetaFYU_imip1_i .* FXV_imip1_i;	% C = bc
    D = D + ThetaFYU_imip1_i .* FYV_imip1_i;	% D = bd
end


%   SVD
[U11,U12, U21, U22, V11, V12, V21, V22] = SVD2D(A, B, C, D);


%   L = UV'
L11 = U11 .* V11 + U12 .* V12;
L12 = U11 .* V21 + U12 .* V22;
L21 = U21 .* V11 + U22 .* V12;
L22 = U21 .* V21 + U22 .* V22;


%   LV
FLVX_1 = L11 .* FXV(:, 1) + L12 .* FYV(:, 1);
FLVY_1 = L21 .* FXV(:, 1) + L22 .* FYV(:, 1);
FLVX_2 = L11 .* FXV(:, 2) + L12 .* FYV(:, 2);
FLVY_2 = L21 .* FXV(:, 2) + L22 .* FYV(:, 2);
FLVX_3 = L11 .* FXV(:, 3) + L12 .* FYV(:, 3);
FLVY_3 = L21 .* FXV(:, 3) + L22 .* FYV(:, 3);
FLVX = [FLVX_1 FLVX_2 FLVX_3];
FLVY = [FLVY_1 FLVY_2 FLVY_3];


%   cot(Theta)LV
FThetaLVX_1 = FCotTheta(:, 1) .* (FLVX_2 - FLVX_3);
FThetaLVY_1 = FCotTheta(:, 1) .* (FLVY_2 - FLVY_3);
FThetaLVX_2 = FCotTheta(:, 2) .* (FLVX_3 - FLVX_1);
FThetaLVY_2 = FCotTheta(:, 2) .* (FLVY_3 - FLVY_1);
FThetaLVX_3 = FCotTheta(:, 3) .* (FLVX_1 - FLVX_2);
FThetaLVY_3 = FCotTheta(:, 3) .* (FLVY_1 - FLVY_2);
MThetaLVX = sparse(Face, Face(:, [2 3 1]), [FThetaLVX_3 FThetaLVX_1 FThetaLVX_2], nV, nV);
MThetaLVY = sparse(Face, Face(:, [2 3 1]), [FThetaLVY_3 FThetaLVY_1 FThetaLVY_2], nV, nV);


end