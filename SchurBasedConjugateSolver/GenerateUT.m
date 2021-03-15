function [Vertex, Face] = GenerateUT(nV)
%
%   this function genertate unstructured triangulation with size nV (using the delaunay method)
%
%   INPUT:  nV - mesh size
%
%   OUTPUT: Vertex - mesh vertices
%           Face - mesh faces


P = gallery('uniformdata',[nV 2],0);
DT = delaunayTriangulation(P);


Vertex = DT.Points;
Vertex(:, 3) = 0;
Face = DT.ConnectivityList;


end