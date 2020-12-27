%%  solving sparse system MX = C directly
%%
%%  Input:
%%          M:	M
%%          C:	C
%%
%%  Output:
%%          X:	X
%%          t:  running time

function [X, t] = DirectSolver(M, C)

%%
tic;
%
X = M \ C;
%
t = toc;

%%
end