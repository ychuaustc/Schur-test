%%  generate 2D rectangular mesh of # nv in 3D space (structured)
%%
%%  Input:  nv: mesh size (constant)
%%  Output: Vertex: vertice matrix (sparse matrix of size vs * 3, vs is mesh size, each row is the x,y,z coordinate
%%          of a vertex in the mesh)
%%          Face: face matrix (sparse matrix of size fs * 3, fs is face size, each row is the vertex index of one mesh
%%          ordered conterclockwisely)

function [Vertex, Face] = GenerateSR(nv)

%%
fprintf('mesh generating...\n');

%%
nx = sqrt(nv);
ny = sqrt(nv);
x = linspace(0, 2, nx);
y = linspace(0, 1, ny);

Vertex = zeros(nv, 3);
Face = zeros(nf, 3);

%%
for i = 1:ny
    Vertex(nx*(i - 1) + 1:nx*i, 1) = x';
    Vertex(nx*(i - 1) + 1:nx*i, 2) = y(i);
end
Vertex(:,3) = 0;

%%
for i = 1:ny - 1
    Face((nx - 1) * (i - 1) + 1:(nx - 1) * i, 1) = nx * (i - 1) + 1:nx * i - 1;
    Face((nx - 1) * (i - 1) + 1:(nx - 1) * i, 2) = nx * (i - 1) + 2:nx * i;
    Face((nx - 1) * (i - 1) + 1:(nx - 1) * i, 3) = nx * i + 2:nx * (i + 1);
    Face((nx - 1) * (i - 1) + 1:(nx - 1) * i, 4) = nx * i + 1:nx * (i + 1) -1;
end

%%
fprintf('mesh generation completed\n\n');
end