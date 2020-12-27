%%  return coefficient matrix for 1D poisson equation (-laplace u + u = f) with p1 element on structured 
%%  triangular mesh ([0, 1] * [0, 1])
%%
%%	Input:
%%          I / B:	inner part / boundary
%%          nI:     size of inner part
%%          nV:     mesh size
%%
%%	Output:
%%          M:      coefficient matrix

function M = Poisson1D_P1_M(I, nI, B, nV)

%%
n = sqrt(nV);
%
lowerB = [2:1:n-1]';
upperB = [(n - 1) * n + 2:1:nV - 1]';
leftB = [n + 1:n:(n - 2) * n + 1]';
rightB = [2 * n:n:(n - 1) * n]';
nlowerB = size(lowerB, 1);
nupperB = size(upperB, 1);
nleftB = size(leftB, 1);
nrightB = size(rightB, 1);
%
leftlowerCo = 1;
rightlowerCo = n;
leftupperCo = (n - 1) * n + 1;
rightupperCo = nV;
Co = [leftlowerCo rightlowerCo leftupperCo rightupperCo]';
nCo = size(Co, 1);
%
BTemp = setdiff(B, Co);
nBTemp = size(BTemp, 1);

%%  -laplace u
%  1 i i
V1_1 = 4 * ones(nI, 1);
V1_2 = 2 * ones(nBTemp, 1);
V1_3 = 1 * ones(nCo, 1);
M1_1 = sparse(I, I, V1_1, nV, nV);
M1_2 = sparse(BTemp, BTemp, V1_2, nV, nV);
M1_3 = sparse(Co, Co, V1_3, nV, nV);
M1 = M1_1 + M1_2 + M1_3;
%	2 i i + 1
V2_1 = (-1) * ones(nI, 1);
V2_2 = -(1 / 2) * ones(nupperB + nlowerB, 1);
V2_3 = (-1) * ones(nleftB, 1);
M2_1 = sparse(I, I + 1, V2_1, nV, nV);
M2_2 = sparse([upperB lowerB], [upperB lowerB] + 1, V2_2, nV, nV);
M2_3 = sparse(leftB, leftB + 1, V2_3, nV, nV);
M2 = M2_1 + M2_2 + M2_3;
M2(leftlowerCo, leftlowerCo + 1) = -(1 / 2);
M2(leftupperCo, leftupperCo + 1) = -(1 / 2);
%	3 i i + n
V3_1 = (-1) * ones(nI, 1);
V3_2 = -(1 / 2) * ones(nleftB + nrightB, 1);
V3_3 = (-1) * ones(nlowerB, 1);
M3_1 = sparse(I, I + n, V3_1, nV, nV);
M3_2 = sparse([leftB rightB], [leftB rightB] + n, V3_2, nV, nV);
M3_3 = sparse(lowerB, lowerB + n, V3_3, nV, nV);
M3 = M3_1 + M3_2 + M3_3;
M3(leftlowerCo, leftlowerCo + n) = -(1 / 2);
M3(rightlowerCo, rightlowerCo + n) = -(1 / 2);

%%  u
%  4 i i
V4_1 = (1 / 2) * ones(nI, 1);
V4_2 = (1 / 4) * ones(nBTemp, 1);
M4_1 = sparse(I, I, V4_1, nV, nV);
M4_2 = sparse(BTemp, BTemp, V4_2, nV, nV);
M4 = M4_1 + M4_2;
M4(leftlowerCo, leftlowerCo) = 1 / 6;
M4(rightlowerCo, rightlowerCo) = 1 / 12;
M4(leftupperCo, leftupperCo) = 1 / 12;
M4(rightupperCo, rightupperCo) = 1 / 6;
%	5 i i + 1
V5_1 = (1 / 12) * ones(nI, 1);
V5_2 = (1 / 24) * ones(nlowerB + nupperB, 1);
V5_3 = (1 / 12) * ones(nleftB, 1);
M5_1 = sparse(I, I + 1, V5_1, nV, nV);
M5_2 = sparse([lowerB upperB], [lowerB upperB] + 1, V5_2, nV, nV);
M5_3 = sparse(leftB, leftB + 1, V5_3, nV, nV);
M5 = M5_1 + M5_2 + M5_3;
M5(leftlowerCo, leftlowerCo + 1) = 1 / 24;
M5(leftupperCo, leftupperCo + 1) = 1 / 24;
%	6 i i + n
V6_1 = (1 / 12) * ones(nI, 1);
V6_2 = (1 / 24) * ones(nleftB + nrightB, 1);
V6_3 = (1 / 12) * ones(nlowerB, 1);
M6_1 = sparse(I, I + n, V6_1, nV, nV);
M6_2 = sparse([leftB rightB], [leftB rightB] + n, V6_2, nV, nV);
M6_3 = sparse(lowerB, lowerB + n, V6_3, nV, nV);
M6 = M6_1 + M6_2 + M6_3;
M6(leftlowerCo, leftlowerCo + n) = 1 / 24;
M6(rightlowerCo, rightlowerCo + n) = 1 / 24;
%	7 i i + n + 1
V7_1 = (1 / 12) * ones(nI, 1);
V7_2 = (1 / 12) * ones(nlowerB + nleftB, 1);
M7_1 = sparse(I, I + n + 1, V7_1, nV, nV);
M7_2 = sparse([lowerB leftB], [lowerB leftB] + n + 1, V7_2, nV, nV);
M7 = M7_1 + M7_2;
M7(leftlowerCo, leftlowerCo + n + 1) = 1 / 12;

%%
MTemp1 = M1 + M4;
MTemp2 = M2 + M3 + M5 + M6 + M7;
M = MTemp1 + MTemp2 + MTemp2';

%%
end