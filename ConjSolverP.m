%%  solve the linear system using conjugate gradient method
%%
%%  Input:  //
%%  Output: //

function [XS, XW, XE, iteration, t] = ConjSolverP(debInd, dewbInd, MS, MSE, MW, MWE, ME, MB, MSB, MWB, MWWB, CS, CW, CE, numDecompose)

%%
tic;

%% compute XE
%
b = CE;
for i = 1:numDecompose
    b = b - MSE{i}' * (MS{i} \ CS{i});
end
if sum(MW) ~= 0
    b = b - MWE' * (MW \ CW);
end
% b = Preconditioner(debInd, dewbInd, MS, MW, MB, MSB, MWB, MWWB, numDecompose, b);
b = PreconditionerW(MW, ME, MWE, b);
% b = PreconditionerW1(debInd, MS, MW, MB, MSB, ME, MWE, numDecompose, b);
% b = PreconditionerW2(MS, MSE, MW, ME, MWE, numDecompose, b);

l = length(b);
XE = zeros(l, 1);
Sx = SchurMultiply(MS, MSE, MW, MWE, ME, XE, numDecompose);
% Sx = Preconditioner(debInd, dewbInd, MS, MW, MB, MSB, MWB, MWWB, numDecompose, Sx);
Sx = PreconditionerW(MW, ME, MWE, Sx);
% Sx = PreconditionerW1(debInd, MS, MW, MB, MSB, ME, MWE, numDecompose, Sx);
% Sx = PreconditionerW2(MS, MSE, MW, ME, MWE, numDecompose, Sx);
r = b - Sx;
p = r;
rsold = r' * r;
%
iteration = 0;
for i = 1:l
    Sp = SchurMultiply(MS, MSE, MW, MWE, ME, p, numDecompose);
%     Sp = Preconditioner(debInd, dewbInd, MS, MW, MB, MSB, MWB, MWWB, numDecompose, Sp);
    Sp = PreconditionerW(MW, ME, MWE, Sp);
%     Sp = PreconditionerW1(debInd, MS, MW, MB, MSB, ME, MWE, numDecompose, Sp);
%     Sp = PreconditionerW2(MS, MSE, MW, ME, MWE, numDecompose, Sp);
    alpha = rsold / (p' * Sp);
    XE = XE + alpha * p;
    r = r - alpha * Sp;
    rsnew = r' * r;
    iteration = iteration + 1;
    if sqrt(rsnew) < 1e-3
        break;
    end
    p = r + (rsnew / rsold) * p;
    rsold = rsnew;
end

%% compute XS with XI
for i = 1:numDecompose
    XS{i} = MS{i} \ (CS{i} - MSE{i} * XE);
end
if sum(MW) ~= 0
    XW = MW \ (CW - MWE * XE);
else
    XW = 0;
end

%%
t = toc;

end