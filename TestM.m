%%  
%%
%%  Input:
%%  Output:        
%%          

function [TM, TMEE_W] = TestM(DE, nV)

%%
TM1 = zeros(nV, nV);
TM2 = zeros(nV, nV);
deInd = find(DE);
notindeInd = setdiff(1:nV, deInd);
n = sqrt(nV);

%%
%   ii
for i = 1:n - 2
    for j = 2:n - 1
        TM1(i * n + j, i * n + j) = 4;
    end
end
for i = 2:n - 1
    TM1(i, i) = 2;
end
for i = (n - 1) * n + 2:n * n - 1
    TM1(i, i) = 2;
end
for i = n + 1:n:(n - 2) * n + 1
    TM1(i, i) = 2;
end
for i = 2 * n:n:(n - 1) * n
    TM1(i, i) = 2;
end
TM1(1, 1) = 1;
TM1(n, n) = 1;
TM1((n - 1) * n + 1, (n - 1) * n + 1) = 1;
TM1(n * n, n * n) = 1;
%   ij
for i = 1:nV
    for j = 1:nV
        if abs(i - j) == n + 1 || abs(i - j) == n - 1
            TM1(i, j) = 0;
        end
        if abs(i - j) == 1 || abs(i - j) == n
            TM1(i, j) = -1;
        end
    end
end
for i = n:n:(n - 1) * n
    TM1(i, i + 1) = 0;
    TM1(i + 1, i) = 0;
    TM1(i, i + n + 1) = 0;
    TM1(i + n + 1, i) = 0;
end
for i = 1:n:(n - 2) * n + 1
    TM1(i, i + n - 1) = 0;
    TM1(i + n - 1, i) = 0;
end
for i = 1:n - 1
    TM1(i, i + 1) = -1 / 2;
    TM1(i + 1, i) = -1 / 2;
end
for i = (n - 1) * n + 1:n * n - 1
    TM1(i, i + 1) = -1 / 2;
    TM1(i + 1, i) = -1 / 2;
end
for i = 1:n:(n - 2) * n + 1
    TM1(i, i + n) = -1 / 2;
    TM1(i + n, i) = -1 / 2;
end
for i = n:n:(n - 1) * n
    TM1(i, i + n) = -1 / 2;
    TM1(i + n, i) = -1 / 2;
end

%%
TMEE_WTemp1 = TM1;
TMEE_WTemp1(notindeInd, :) = 0;
TMEE_WTemp1(:, notindeInd) = 0;
%   ii
for i = n + 3:n:(n - 2) * n + 3
    TMEE_WTemp1(i, i) = 2;
end
for i = (n + 1) / 2:n:(n - 1) * n + (n + 1) / 2
    TMEE_WTemp1(i, i) = 0;
end
for i = n + n - 2:n:(n - 2) * n + n - 2
    TMEE_WTemp1(i, i) = 2;
end
TMEE_WTemp1(3, 3) = 1;
TMEE_WTemp1(n - 2, n - 2) = 1;
TMEE_WTemp1((n - 1) * n + 3, (n - 1) * n + 3) = 1;
TMEE_WTemp1(n * n - 2, n * n - 2) = 1;
%   ij
for i = 3:n:(n - 2) * n + 3
    TMEE_WTemp1(i, i + 1) = -1 / 2;
    TMEE_WTemp1(i + 1, i) = -1 / 2;
end
for i = (n + 1) / 2:n:(n - 2) * n + (n + 1) / 2
    TMEE_WTemp1(i, i + 1) = 0;
    TMEE_WTemp1(i + 1, i) = 0;
end
for i = n - 2:n:(n - 2) * n + n - 2
    TMEE_WTemp1(i, i + 1) = -1 / 2;
    TMEE_WTemp1(i + 1, i) = -1 / 2;
end
%
TMEE_W1 = TMEE_WTemp1(deInd, deInd);

%%
TM1 = 2 * TM1;
TMEE_W1 = 2 * TMEE_W1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%   ii
for i = 1:2:nV
    TM2(i, i) = 2 / 3;
end
for i = 2:2:nV - 1
    TM2(i, i) = 1 / 3;
end
for i = 2:2:n - 1
    TM2(i, i) = 1 / 6;
end
for i = (n - 1) * n + 2:2:nV - 1
    TM2(i, i) = 1 / 6;
end
for i = n + 1:2 * n:(n - 2) * n + 1
    TM2(i, i) = 1 / 6;
end
for i = 2 * n:2 * n:(n - 1) * n
    TM2(i, i) = 1 / 6;
end
for i = 3:2:n - 2
    TM2(i, i) = 1 / 3;
end
for i = (n - 1) * n + 3:2:nV - 2
    TM2(i, i) = 1 / 3;
end
for i = 2 * n + 1:2 * n:(n - 3) * n + 1
    TM2(i, i) = 1 / 3;
end
for i = 3 * n:2 * n:(n - 2) * n
    TM2(i, i) = 1 / 3;
end
TM2(1, 1) = 1 / 6;
TM2(n, n) = 1 / 6;
TM2((n - 1) * n + 1, (n - 1) * n + 1) = 1 / 6;
TM2(n * n, n * n) = 1 / 6;
%   ij
for i = 1:nV
    for j = 1:nV
        if abs(i - j) == 1 || abs(i - j) == n
            TM2(i, j) = 1 / 12;
        end
    end
end
for i = 1:2:(n - 1) * n -1
    TM2(i, i + n + 1) = 1 / 12;
end
for i = 3:2:(n - 1) * n -1
    TM2(i, i + n - 1) = 1 / 12;
end
for i = n:n:(n - 1) * n
    TM2(i, i + 1) = 0;
    TM2(i + 1, i) = 0;
    TM2(i, i + n + 1) = 0;
    TM2(i + n + 1, i) = 0;
end
for i = 1:n:(n - 2) * n + 1
    TM2(i, i + n - 1) = 0;
    TM2(i + n - 1, i) = 0;
end
for i = 1:n - 1
    TM2(i, i + 1) = 1 / 24;
    TM2(i + 1, i) = 1 / 24;
end
for i = (n - 1) * n + 1:n * n - 1
    TM2(i, i + 1) = 1 / 24;
    TM2(i + 1, i) = 1 / 24;
end
for i = 1:n:(n - 2) * n + 1
    TM2(i, i + n) = 1 / 24;
    TM2(i + n, i) = 1 / 24;
end
for i = n:n:(n - 1) * n
    TM2(i, i + n) = 1 / 24;
    TM2(i + n, i) = 1 / 24;
end

%%
TMEE_WTemp2 = TM2;
TMEE_WTemp2(notindeInd, :) = 0;
TMEE_WTemp2(:, notindeInd) = 0;
%   ii
for i = n + 3:2 * n:(n - 2) * n + 3
    TMEE_WTemp2(i, i) = 1 / 6;
end
for i = 2 * n + 3:2 * n:(n - 3) * n + 3
    TMEE_WTemp2(i, i) = 1 / 3;
end
for i = (n + 1) / 2:n:(n - 1) * n + (n + 1) / 2
    TMEE_WTemp2(i, i) = 0;
end
for i = 2 * n - 2:2 * n:(n - 1) * n - 2
    TMEE_WTemp2(i, i) = 1 / 6;
end
for i = 3 * n - 2:2 * n:(n - 2) * n - 2
    TMEE_WTemp2(i, i) = 1 / 3;
end
TMEE_WTemp2(3, 3) = 1 / 6;
TMEE_WTemp2(n - 2, n - 2) = 1 / 6;
TMEE_WTemp2((n - 1) * n + 3, (n - 1) * n + 3) = 1 / 6;
TMEE_WTemp2(n * n - 2, n * n - 2) = 1 / 6;
%   ij
for i = 3:n:(n - 2) * n + 3
    TMEE_WTemp2(i, i + 1) = 1 / 24;
    TMEE_WTemp2(i + 1, i) = 1 / 24;
end
for i = (n + 1) / 2:n:(n - 2) * n + (n + 1) / 2
    TMEE_WTemp2(i, i + 1) = 0;
    TMEE_WTemp2(i + 1, i) = 0;
end
for i = n - 2:n:(n - 2) * n + n - 2
    TMEE_WTemp2(i, i + 1) = 1 / 24;
    TMEE_WTemp2(i + 1, i) = 1 / 24;
end
%
TMEE_W2 = TMEE_WTemp2(deInd, deInd);

%%
TM2 = 2 * TM2;
TMEE_W2 = 2 * TMEE_W2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TM = TM1 + TM2;
TMEE_W = TMEE_W1 + TMEE_W2;

%%
end