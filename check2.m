%%  check the validation of the wirebasket set
%%
%%  Input:  nV: mesh size (const)
%%          DS: list of subdomain for the decomposition (tuple of sparse nV * 1 matrix)
%%          B: boundary (column vector)
%%          DE: edges in the decomposition (sparse nV * 1 matrix)
%%          DW: vertice (wirebasket) for the decomposition (sparse nV * 1 matrix)
%%          numDecompose: number of decomposition (constant)
%%  Output: che2: check result (constant)

function [che2] = check2(nV, DS, B, DE, DW, numDecompose)

%%

a = sum(double(DE < 0));
b = sum(double(DW < 0));
c = sum(double(DE & DW));
d = DE;
d = d + DW;
e = ones(nV, 1);
e(B) = 0;
for i = 1:numDecompose
    e = e - DS{i};
end
f = d - e;
g = sum(f ~= 0);

h = [a, b, c, g];

che2 = sum(h ~= 0);

end