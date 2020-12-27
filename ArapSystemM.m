%%  return M for MX = C
%%
%%  Input:
%%          MCotTheta:  
%%          nV:         
%%          epsLSystem: 
%%  Output:
%%          M:          
%%          

function [M] = ArapSystemM(MCotTheta, nV)

%%
%
M = sparse(1:nV, 1:nV, sum(MCotTheta + MCotTheta', 2), nV, nV) - (MCotTheta + MCotTheta');

%%
end