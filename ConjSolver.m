%%  solve the linear system using conjugate gradient method
%%
%%  Input:  //
%%  Output: //

function [XS, XW, XE, iteration, t] = ConjSolver(MS, MSE, MW, MWE, ME, CS, CW, CE, numDecompose)

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

l = length(b);
XE = zeros(l, 1);
Sx = SchurMultiply(MS, MSE, MW, MWE, ME, XE, numDecompose);
r = b - Sx;
p = r;
rsold = r' * r;
%
iteration = 0;
for i = 1:l
    Sp = SchurMultiply(MS, MSE, MW, MWE, ME, p, numDecompose);
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