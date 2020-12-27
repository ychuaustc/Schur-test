%%  check the validation of the wirebasket set
%%
%%  Input:
%%          nV:             mesh size
%%          dsall:          all subdomains for the decomposition
%%          DE:             edges in the decomposition
%%          DW:             wirebasket set in the decomposition
%%          B:              boundary
%%          numDecompose:   number of decompositions
%%  Output:

function [] = DecomposeUTCheck2(nV, dsall, DE, DW)

%%

a = sum(double(DE < 0));
b = sum(double(DW < 0));
c = sum(double(DE & DW));
d = DE + DW;
e = ones(nV, 1);
e(dsall) = 0;
f = d - e;
g = sum(f);
%
h = [a, b, c, g];
%
che2 = sum(h);

%%
if che2 == 0
    fprintf('wirebasket set computed\n\n');
else
    fprintf('wirebasket computation failed\n');
    quit(1);
end

%%
end