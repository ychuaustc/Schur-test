function [] = PlotMesh(Vertex, Face, k)
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
switch k
    case 2
        xlabel(['Direct Solver'],'visible','on');
    case 3
        xlabel(['Conjugate Gradient Solver'],'visible','on');
    case 4
        xlabel(['Preconditioned Conjugate Gradient Solver'],'visible','on');
    otherwise
        quit(1);
end

hold on;

set(trimesh(Face, Vertex(:, 1), Vertex(:, 2)), 'color', 'k');

hold off;

box off;


end