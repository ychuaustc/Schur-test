%%  return U = PM^-1 * V, where PM is the preconditioner for Schur matrix
%%
%%  Input:
%%  Output:

function U = PreconditionerW(MW, ME, MWE, V)

%%
if sum(MW) ~= 0
    U = (ME - MWE' * (MW \ MWE)) \ V;
else
    U = V;
end

end