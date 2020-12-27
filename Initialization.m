%%  return a initial parameterization
%%
%%  Input:
%%          MC:         connection matrix
%%          I / B:      inner part / boundary
%%          nB:         size of boundary
%%          nV:         mesh size
%%  Output:
%%          Vertex:     parameterized vertice matrix

function [Vertex] = Initialization(MC, I, B, nB, nV)

%%
Vertex = zeros(nV, 2);

%%   fix the boundary
%
BX = cos([0:2 * pi / nB:2 * pi])';
BY = sin([0:2 * pi / nB:2 * pi])';
BX = BX(1:1:nB);
BY = BY(1:1:nB);
%
Vertex(B, 1) = BX;
Vertex(B, 2) = BY;

%%  fix the inside part
%
M = sparse(1:nV, 1:nV, sum(MC, 2), nV, nV) - MC;
MI = M(I, I);
MB = MC(I, B);
%
IX = MI \ (MB * BX);
IY = MI \ (MB * BY);
%
Vertex(I, 1) = IX;
Vertex(I, 2) = IY;

%%
end