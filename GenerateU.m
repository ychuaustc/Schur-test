%%  generate 2D triangular mesh of # nv in 3D space (unstructed)
%%
%%  Input:  nv: mesh size (constant)
%%  Output: Vertex: vertice matrix (sparse matrix of size vs * 3, vs is mesh size, each row is the x,y,z coordinate
%%          of a vertex in the mesh)
%%          Face: face matrix (sparse matrix of size fs * 3, fs is face size, each row is the vertex index of one mesh
%%          ordered conterclockwisely)

function [Vertex, Face] = GenerateU(nv)

%%
fprintf('mesh generating...\n');

%%
P = gallery('uniformdata',[nv 2],0);
DT = delaunayTriangulation(P);

%%
Vertex = DT.Points;
Vertex(:, 3) = 0;
Face = DT.ConnectivityList;

%%
fprintf('mesh generation completed\n\n');
end