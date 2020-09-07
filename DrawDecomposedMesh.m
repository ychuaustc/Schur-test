%%  draw mesh
%%  Input:  Vertex: vertice matrix (sparse matrix, vertice size * 3)
%%          Face: face matrix (sparse matrix, face size * 3)
%%          DV: list of vertice in the decomposition (tuple of sparse nV * 1 matrix)
%%          nV: mesh size (const)
%%          numDecompose: number of decomposition (constant)
%%  Output: //

function [] = DrawDecomposedMesh(Vertex, Face, DV, nV, numDecompose)
    
%%
figure; set(gcf, 'Units', 'normalized', 'Position', [0.05, 0.05, .8, .8]);

%%
plot([0, 1], [0, 0], 'k');hold on; plot([0, 1], [1, 1], 'k');hold on;
plot([0, 0], [0, 1], 'k');hold on; plot([1, 1], [0, 1], 'k');hold on;
axis off; axis equal; title('Decomposed mesh');

%%
for i = 1:numDecompose
    facecolor = rand(1, 3);
    tempVertice = find(ones(nV, 1) - DV{i})';
    tempFace = Face;
    for j = 1:size(tempVertice, 2)
        for k = 1:3
            tempFace(find(tempFace(:, k) == tempVertice(j))', :) = [];
        end
    end
    
    trimesh(tempFace, Vertex(:, 1), Vertex(:, 2), Vertex(:, 3), 'facecolor', facecolor);
    hold on;
end

end