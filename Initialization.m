function [Vertex, X, Y] = Initialization(MC, I, B, nB, nV)
%
%   this function computes a intial prameterization using the Tutte method
%
%   INPUT:  MC - adjacent matrix
%           I - inner part
%           B - boundary
%           nB - boundary size
%           nV - mesh size
%
%   OUTPUT: Vertex - mesh vertices
%           X - the x-coordinates
%           Y - the y-coordinates


Vertex = zeros(nV, 2);


%   fix the boundary
BX = cos([0:2 * pi / nB:2 * pi])';
BY = sin([0:2 * pi / nB:2 * pi])';
BX = BX(1:1:nB);
BY = BY(1:1:nB);
Vertex(B, 1) = BX;
Vertex(B, 2) = BY;


M = sparse(1:nV, 1:nV, sum(MC, 2), nV, nV) - MC;
MI = M(I, I);
MB = MC(I, B);


IX = MI \ (MB * BX);
IY = MI \ (MB * BY);


Vertex(I, 1) = IX;
Vertex(I, 2) = IY;


X = Vertex(:, 1);
Y = Vertex(:, 2);


end