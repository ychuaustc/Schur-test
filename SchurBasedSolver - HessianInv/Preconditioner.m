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


U = (MEE_W - MWE' * (MWW \ MWE)) \ V;


end