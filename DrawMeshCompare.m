%%  draw mesh
%%  Input:  Vertex1: vertice 1 (sparse matrix, vertice size * 3)
%%          Vertex2: vertice 2 (sparse matrix, vertice size * 3)
%%          Vertex3: vertice 3 (sparse matrix, vertice size * 3)
%%          Vertex4: vertice 4 (sparse matrix, vertice size * 3)
%%          Face: face (sparse matrix, face size * 3)
%%  Output: draw original and parameterized mesh on [0, 1] * [0, 1]

function [] = DrawMeshCompare(Vertex1, Vertex2, Vertex3, Vertex4, Face)
    
%%
figure; set(gcf, 'Units', 'normalized', 'Position', [0.05, 0.05, .8, .8]);

%%
subplot(221);
plot([0, 1], [0, 0], 'k');hold on; plot([0, 1], [1, 1], 'k');hold on;
plot([0, 0], [0, 1], 'k');hold on; plot([1, 1], [0, 1], 'k');hold on;
trimesh(Face, Vertex1(:, 1), Vertex1(:, 2), Vertex1(:, 3), 'edgecolor', 'k');
axis off; axis equal; title('Parameterization');

%%
subplot(222);
plot([0, 1], [0, 0], 'k');hold on; plot([0, 1], [1, 1], 'k');hold on;
plot([0, 0], [0, 1], 'k');hold on; plot([1, 1], [0, 1], 'k');hold on;
trimesh(Face, Vertex2(:, 1), Vertex2(:, 2), Vertex2(:, 3), 'edgecolor', 'k');
axis off; axis equal; title('Parameterization using direct solver');

%%
subplot(223);
plot([0, 1], [0, 0], 'k');hold on; plot([0, 1], [1, 1], 'k');hold on;
plot([0, 0], [0, 1], 'k');hold on; plot([1, 1], [0, 1], 'k');hold on;
trimesh(Face, Vertex3(:, 1), Vertex3(:, 2), Vertex3(:, 3), 'edgecolor', 'k');
axis off; axis equal; title('Parameterization using conjugate solver');

%%
subplot(224);
plot([0, 1], [0, 0], 'k');hold on; plot([0, 1], [1, 1], 'k');hold on;
plot([0, 0], [0, 1], 'k');hold on; plot([1, 1], [0, 1], 'k');hold on;
trimesh(Face, Vertex4(:, 1), Vertex4(:, 2), Vertex4(:, 3), 'edgecolor', 'k');
axis off; axis equal; title('Parameterization using conjugate solver with preconditioner');

end