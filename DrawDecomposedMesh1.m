%%  draw mesh
%%  Input:
%%          Vertex:         vertice matrix
%%          Face:           face matrix
%%          DV:             list of vertice in the decomposition
%%          nV:             mesh size
%%          numDecompose:	number of decomposition
%%  Output: 

function [] = DrawDecomposedMesh1(Vertex, Face, DV, nV, numDecompose)
    
%%
subplot(2, 2, 1);

%%
plot([0, 1], [0, 0], 'k');hold on; plot([0, 1], [1, 1], 'k');hold on;
plot([0, 0], [0, 1], 'k');hold on; plot([1, 1], [0, 1], 'k');hold on;
axis off; axis equal; title('Decomposed mesh');

%%
hold on;

set(trimesh(Face, Vertex(:, 1), Vertex(:, 2)), 'color', 'k');

hold off;

Z = zeros(nV, 1);
for i = 1:numDecompose
% for i = 1:numDecompose + 1
    facecolor = rand(1, 3);
    tempVertice = find(ones(nV, 1) - DV{i})';
    tempFace = Face;
    for j = 1:size(tempVertice, 2)
        for k = 1:3
            tempFace(find(tempFace(:, k) == tempVertice(j))', :) = [];
        end
    end
    %
    trimesh(tempFace, Vertex(:, 1), Vertex(:, 2), Z, 'facecolor', facecolor);
    hold on;
end

%%
end