%%  return U = Schur matrix * V
%%
%%  Input:  MS
%%          MSE
%%          MW
%%          MWE
%%          ME
%%          V
%%          numDecompose: number of decomposition (constant)
%%  Output: U: Schur matrix * V (column vector)

function U = SchurMultiply(MS, MSE, MW, MWE, ME, V, numDecompose)

%%
U = ME * V;
for i = 1:numDecompose
    U = U - MSE{i}' * (MS{i} \ (MSE{i} * V));
end
if sum(MW) ~= 0
    U = U - MWE' * (MW \ (MWE * V));
end

end