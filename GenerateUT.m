%%  generate 2D triangular mesh of # nV (unstructed)
%%
%%  Input:
%%          nV:     mesh size
%%  Output:
%%          Vertex: vertice matrix
%%          Face:   face matrix

function [Vertex, Face] = GenerateUT(nV)

%%
fprintf('mesh generating...\n');

%%
P = gallery('uniformdata',[nV 2],0);
DT = delaunayTriangulation(P);

%%
Vertex = DT.Points;
Face = DT.ConnectivityList;

%%
fprintf('mesh generation completed\n\n');

%%
end