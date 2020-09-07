%%  return U = PM^-1 * V, where PM is the preconditioner for Schur matrix
%%
%%  Input:
%%  Output:

function U = PreconditionerW2(MS, MSE, MW, ME, MWE, numDecompose, V)

%%
l = size(V, 1);
U = zeros(l, 1);

%%
% if sum(MW) ~= 0
%     w = 1 / (numDecompose + 1);
% else
%     w = 1/ numDecompose;
% end
w = 1;

%%
for i = 1:numDecompose
    U = U + w * w * (ME - MSE{i}' * (MS{i} \ MSE{i})) \ V;
end

%%
if sum(MW) ~= 0
    U = U + w * w * (ME - MWE' * (MW \ MWE)) \ V;
end

end