%%  return energy of U and LX
%%
%%  Input:
%%          FCotTheta:  
%%          FLVX:       
%%          FLVY:       
%%          FXU:        
%%          FYU:        
%%  Output:
%%          EnUL:       

function [EnUL] = EnergyUL(FCotTheta, FLVX, FLVY, FXU, FYU)

%%
%
EnUL = 0;
%
for i = 1:3
    ip1 = mod(i, 3) + 1;    % (i + 1) mod 3
    %
    CXi = (FXU(:, i) - FXU(:, ip1)) - (FLVX(:, i) - FLVX(:, ip1));
    CYi = (FYU(:, i) - FYU(:, ip1)) - (FLVY(:, i) - FLVY(:, ip1));
    Ei = 0.5 * sum(FCotTheta(:, i) .* (CXi .* CXi + CYi .* CYi));
    EnUL = EnUL + Ei;
end

%%
end