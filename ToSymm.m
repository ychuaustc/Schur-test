%%  turn a non symmetric matrix to symmtric
%%
%%  Input:
%%          
%%  Output:
%%          

function [MSymm] = ToSymm(M, nV)

%%
MTemp = M';
[x, y] = find(M - MTemp);
z = ((M - MTemp) ~= 0);
Mdiff = MTemp(z);
MDiff = sparse(x, y, Mdiff, nV,nV);
MSymm = M + MDiff;

%%
end