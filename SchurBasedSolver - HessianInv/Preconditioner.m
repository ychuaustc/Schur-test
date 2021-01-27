function U = Preconditioner(MWW, MWE, MEE_W, V)
%
%   this function applies the Neumann-Dirichlet preconditioner for the Schur system (M^(-1)V)
%
%   INPUT:  MWW - the wirebasket matrix for the Schur system
%           MWE - the wirebasket-edge matrix for the Schur system
%           MEE_W - the wirebasket matrix block
%           V - the vector to be preconditioned
%
%   OUTPUT: U - the preconditioned vector M^(-1)V

% tic;
% U = (MEE_W - MWE' * (MWW \ MWE)) \ V;
nw = size(MWW, 1);
nb = size(MEE_W, 1);
zw = zeros(nw, 1);
U = [MWW MWE;MWE' MEE_W] \ [zw;V];
U = U(nw + 1:nw + nb);
% t = toc;

end