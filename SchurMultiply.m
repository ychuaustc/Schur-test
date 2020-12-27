%%  return U = Schur matrix S * V
%%
%%  Input:
%%          MSS:            matrix block for subdomain * subdomain (tuple)
%%          MSE:            matrix block for subdomain * edge (tuple)
%%          MWW:            matrix block for wirebasket * wirebasket
%%          MWE:            matrix block for wirebasket * edge
%%          MEE:            matrix block for edge * edge
%%          V:              mutiplier vector
%%          numDecompose:	number of decomposition
%%  Output: U:              Schur matrix S * V

function U = SchurMultiply(MSS, MSE, MWW, MWE, MEE, V, numDecompose)

%%
U = MEE * V;
for i = 1:numDecompose
    U = U - MSE{i}' * (MSS{i} \ (MSE{i} * V));
end
U = U - MWE' * (MWW \ (MWE * V));

end