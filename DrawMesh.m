function [] = DrawMesh(Vertex, Face, nV)
%
%   this function draws the parameterized mesh
%
%   INPUT:  Vertex - mesh vertices
%           Face - mesh faces
%           nV - mesh size
%
%   OUTPUT: 
    

figure; set(gcf, 'Units', 'normalized', 'Position', [0.05, 0.05, .8, .8]);


plot([0, 1], [0, 0], 'k');hold on;plot([0, 1], [1, 1], 'k');hold on;
plot([0, 0], [0, 1], 'k');hold on;plot([1, 1], [0, 1], 'k');hold on;
axis off; axis equal; title('Parameterized mesh');


Z = zeros(nV, 1);
trimesh(Face, Vertex(:, 1), Vertex(:, 2), Z, 'edgecolor', 'k');


end