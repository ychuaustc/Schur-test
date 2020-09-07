%%  generate 2D triangular mesh of # nv in 3D space (structured)
%%
%%  Input:  nv: mesh size (constant)
%%  Output: Vertex: vertice matrix (sparse matrix of size vs * 3, vs is mesh size, each row is the x,y,z coordinate
%%          of a vertex in the mesh)
%%          Face: face matrix (sparse matrix of size fs * 3, fs is face size, each row is the vertex index of one mesh
%%          ordered conterclockwisely)

function [Vertex, Face] = GenerateS(nv)

%%
fprintf('mesh generating...\n');

%%
nx = sqrt(nv);
ny = sqrt(nv);
nfh = (nx - 1) * (ny - 1);
nf = nfh * 2;
x = linspace(0, 2, nx);
y = linspace(0, 1, ny);

Vertex = zeros(nv, 3);
Face1 = zeros(nfh, 3);
Face2 = zeros(nfh, 3);
Face = zeros(nf, 3);

%%
for i = 1:ny
    Vertex(nx*(i - 1) + 1:nx*i, 1) = x';
    Vertex(nx*(i - 1) + 1:nx*i, 2) = y(i);
end
Vertex(:,3) = 0;

%%
for i = 1:ny - 1
    Face1((nx - 1) * (i - 1) + 1:(nx - 1) * i, 1) = nx * (i - 1) + 1:nx * i - 1;
    Face1((nx - 1) * (i - 1) + 1:(nx - 1) * i, 2) = nx * i + 2:nx * (i + 1);
    Face1((nx - 1) * (i - 1) + 1:(nx - 1) * i, 3) = nx * i + 1:nx * (i + 1) -1;
    Face2((nx - 1) * (i - 1) + 1:(nx - 1) * i, 1) = nx * (i - 1) + 1:nx * i - 1;
    Face2((nx - 1) * (i - 1) + 1:(nx - 1) * i, 2) = nx * (i - 1) + 2:nx * i;
    Face2((nx - 1) * (i - 1) + 1:(nx - 1) * i, 3) = nx * i + 2:nx * (i + 1);
end
Face(1:nfh, :) = Face1;
Face(nfh + 1:nf, :) = Face2;

%%
fprintf('mesh generation completed\n\n');
end