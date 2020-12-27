%%  return the information of the mesh
%%
%%  Input:
%%          Vertex:     vertice matrix
%%          Face:       face matrix
%%          nV:         mesh size
%%  Output:
%%          nF:         size of faces
%%          I / B:      inner part / boundary
%%          nI / nB:    size of inner part / boundary
%%          MC:         connection matrix

function [nF, I, nI, B, nB, MC] = MeshInfoUT(Vertex, Face, nV)

%%
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
fprintf('mesh information:\nnumber of vertice: %d \nnumber of face: %d\n\n\n\n', nV, nF);

%%
end