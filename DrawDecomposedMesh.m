function [] = DrawDecomposedMesh(Vertex, Face, DV, nV, nF, numDecompose)
%
%   this function draws all partitions of the decomposed mesh
%
%   INPUT:  Vertex - mesh vertices
%           Face - mesh faces
%           DV - decomposition list
%           nV - mesh size
%           nF - number of faces
%           numDecompose - decomposition number
%
%   OUTPUT: 
    

figure; set(gcf, 'Units', 'normalized', 'Position', [0.05, 0.05, .8, .8]);


plot([0, 1], [0, 0], 'k');hold on; plot([0, 1], [1, 1], 'k');hold on;
plot([0, 0], [0, 1], 'k');hold on; plot([1, 1], [0, 1], 'k');hold on;
axis off; axis equal; title('Decomposed mesh');


Z = zeros(nV, 1);
for i = 1:numDecompose
    facecolor = rand(1, 3);
    tempVertice = find(ones(nV, 1) - DV{i})';
    tempFace = Face;
    for j = 1:size(tempVertice, 2)
        for k = 1:3
            % mark the faces which contains the vertice that is not in the current decomposition
            tempFace(:, k) = tempFace(:, k) .* double(tempFace(:, k) ~= tempVertice(j));
        end
    end
    v = setdiff(1:nF, find((tempFace(:, 1) .* tempFace(:, 2)) .* tempFace(:, 3))');
    tempFace(v, :) = [];    % remove all marked faces
    
    trimesh(tempFace, Vertex(:, 1), Vertex(:, 2), Z, 'facecolor', facecolor);
    hold on;
end


end