%%  return the information of the mesh
%%
%%  Input:  Vertex: vertice matrix (sparse matrix of size vs * 3, vs is mesh size, each row is the x,y,z coordinate
%%          of a vertex in the mesh)
%%          Face: face matrix (sparse matrix of size fs * 3, fs is face size, each row is the vertex index of one mesh
%%          ordered conterclockwisely)
%%  Output: nV / nF: size of vertice / face (costant)
%%          B / I: boundary / inner part of mesh (column vector)
%%          nB / nI: size of B / I
%%          MC: connection matrix (sparse nV * nV matrix)

function [nV, nF, B, nB, I, nI, MC] = MeshInfo(Vertex, Face)

%%
nV = size(Vertex, 1);
nF = size(Face, 1);

%% generate the triangulation meshes
TR = triangulation(Face,Vertex);

%% find boundary vertex & inner vertex
B = freeBoundary(TR);
B = B(:, 1);
nB = size(B, 1);

I = setdiff(1:nV, B)';
nI = size(I, 1);

%% compute connection matrix
MCtemp = sparse(Face, Face(:, [2 3 1]), 1, nV, nV);
MC = double(MCtemp | MCtemp');
fprintf('mesh information:\nnumber of vertis: %d \nnumber of face: %d\n\n\n\n', nV, nF);

end