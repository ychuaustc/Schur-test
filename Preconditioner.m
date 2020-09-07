%%  return U = PM^-1 * V, where PM is the preconditioner for Schur matrix
%%
%%  Input:
%%  Output:

function U = Preconditioner(debInd, dewbInd, MS, MW, MB, MSB, MWB, MWWB, numDecompose, V)

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
    Utemp = zeros(l, 1);
    Utemp(debInd{i}) = (MB{i} - MSB{i}' * (MS{i} \ MSB{i})) \ V(debInd{i});
    U = U + w * w * Utemp;
end

%%
if sum(MW) ~= 0
    Utemp = zeros(l, 1);
    Utemp(dewbInd) = (MWB - MWWB' * (MW \ MWWB)) \ V(dewbInd);
    U = U + w * w * Utemp;
end

end