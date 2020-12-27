%%  generate 2D triangular mesh of # nV (structured)
%%
%%  Input:
%%          nV:       mesh size
%%  Output:
%%          Vertex:   vertice matrix
%%          Face:     face matrix

function [Vertex, Face] = GenerateST(nV)

%%
fprintf('mesh generating...\n');

%%
nx = sqrt(nV);
ny = sqrt(nV);
nfh = (nx - 1) * (ny - 1);
nf = nfh * 2;
x = linspace(0, 1, nx);
y = linspace(0, 1, ny);

%%
Vertex = zeros(nV, 3);
Face = zeros(nf, 3);

%%
for i = 1:ny
    Vertex(nx*(i - 1) + 1:nx*i, 1) = x';
    Vertex(nx*(i - 1) + 1:nx*i, 2) = y(i);
end
Vertex(:,3) = 0;

%%
nx1 = nx - 1;
ny1 = ny - 1;
L = zeros(nx1 * ny1, 1);
for i = 1:1:ny1
    L((i - 1) * nx1 + 1:1:i * nx1) = (i - 1) * nx + 1:1:i * nx - 1;
end
Face(1:2:nf - 1, [1 2 3]) = [L L + nx + 1 L + nx];
Face(2:2:nf, [1 2 3]) = [L L + 1 L + nx + 1];

%%
fprintf('mesh generation completed\n\n');

%%
end