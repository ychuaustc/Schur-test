%%  
%%
%%  Input:
%%          
%%  Output:
%%          

function [DS, DB, DE] = EconfigS(DV, DE, DB, B, nV, nv, numDecompose)
%%
%
DE(B) = 0;
sqnV = sqrt(nV);
sqnv = sqrt(nv);
DE(sqnv:sqnv - 1:sqnV - sqnv + 1) = 1;
DE((sqnV - 1) * sqnV + sqnv:sqnv - 1:nV - sqnv + 1) = 1;
DE((sqnv - 1) * sqnV + 1:(sqnv - 1) * sqnV:nV - sqnv * sqnV + 1) = 1;
DE(sqnv * sqnV:(sqnv - 1) * sqnV:nV - (sqnv - 1) * sqnV) = 1;
%
for i = 1:numDecompose
    DB{i} = DB{i} & DE;
    DS{i} = DV{i} - DB{i};
end

%%
end