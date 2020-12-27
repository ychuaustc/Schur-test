%%	return U = P^-1 * V for given V, where P is the preconditioner for Schur matrix S (temp version, for test)
%%
%%  Input:
%%          MWW:            matrix block for wirebasket * wirebasket
%%          MWE:            matrix block for wirebasket * edge
%%          MEE:            matrix block for edge * edge
%%          V:              vector
%%  Output: U:              P^-1 * V, P is the preconditioner for Schur matrix S

function U = Preconditioner(MWW, MWE, MEE_W, LW, LE, V)

%%
nw = size(MWW, 1);
nb = size(MEE_W, 1);
ze = zeros(nw, 1);
U = [MWW MWE LW;MWE' MEE_W LE;LW' LE' 0] \ [ze;V;0];
U = U(nw + 1:nw + nb);
%
% U = V;
%
% U = (MEE_W - MWE' * (MWW \ MWE)) \ V;

end