function [] = PlotMesh1(Vertex, Face, k)
%
%   this function draws the parameterized mesh
%
%   INPUT:  Vertex - mesh vertices
%           Face - mesh faces
%           nV - mesh size
%           k - plot position
%
%   OUTPUT: 
    

subplot(2, 2, k);

axislim = [min(Vertex);max(Vertex)];
axis(axislim(:)');

axis off;
axis equal;

hold on;

set(trimesh(Face, Vertex(:, 1), Vertex(:, 2)), 'color', 'k');

hold off;

box off;


end