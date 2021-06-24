function [nF, I, nI, B, nB, MC] = MeshInfo(Vertex, Face, nV)
%
%   this function gives the details of the triangulation, such as faces numbers, adjacent matrix and
%   so on
%
%   INPUT:  Vertex - mesh vertices
%           Face - mesh faces
%           nV - mesh size
%
%   OUTPUT: nF - number of faces
%           I - inner part
%           nI - number of inner part
%           B - boundary
%           nB - number of boundary
%           MC - adjacent matrix


nF = size(Face, 1);


TR = triangulation(Face,Vertex);


B = freeBoundary(TR);
B = B(:, 1);
% B = findBoundary(Vertex, Face);
nB = size(B, 1);

I = setdiff(1:nV, B)';
nI = size(I, 1);


MCtemp = sparse(Face, Face(:, [2 3 1]), 1, nV, nV);
MC = double(MCtemp | MCtemp');
fprintf('mesh information:\nthe number of vertices: %d \nthe number of faces: %d\n\n\n', nV, nF);


end