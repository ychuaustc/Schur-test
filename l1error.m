%%  compute the L1-error of two vertice set
%%
%%  Input:  vertex1: vertice matrix 1 (sparse matrix, vertis size * 3)
%%          vertex2: vertice matrix 2 (sparse matrix, vertis size * 3)
%%          nV: mesh size (constant)
%%  Output: l1error: l1-error (constant)

function [l1error] = l1error(vertex1, vertex2, nV)

l1error = sum(sum(abs(vertex1 - vertex2))) / (nV * 3);

end