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
Vertex(:, 3) = Vertex(:, 1).^2 + Vertex(:, 2).^2;
Face = DT.ConnectivityList;


end