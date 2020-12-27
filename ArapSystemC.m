%%  return C for MX = C
%%
%%  Input:
%%          MThetaLV_ij:    
%%          MThetaLV_ji:    
%%  Output:
%%          C:              

function [C] = ArapSystemC(MThetaLV)

%%
C = sum(MThetaLV, 2) - sum(MThetaLV, 1)';

%%
end