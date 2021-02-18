function [Vertex, Face] = GenerateST(nV)
%
%   this function genertate structured triangulation with size nV
%
%   INPUT:  nV - mesh size
%
%   OUTPUT: Vertex - mesh vertices
%           Face - mesh faces


fprintf('mesh generating...\n');


nX = sqrt(nV);
nY = sqrt(nV);
nFhalf = (nX - 1) * (nY - 1);
nF = nFhalf * 2;
x = linspace(0, 1, nX);
y = linspace(0, 1, nY);


Vertex = zeros(nV, 3);
Face = zeros(nF, 3);


for i = 1:nY
    Vertex(nX*(i - 1) + 1:nX*i, 1) = x';
    Vertex(nX*(i - 1) + 1:nX*i, 2) = y(i);
end
Vertex(:, 3) = 0;


nx1 = nX - 1;
ny1 = nY - 1;
L = zeros(nx1 * ny1, 1);
for i = 1:1:ny1
    L((i - 1) * nx1 + 1:1:i * nx1) = (i - 1) * nX + 1:1:i * nX - 1;
end
Face(1:2:nF - 1, [1 2 3]) = [L L + nX + 1 L + nX];
Face(2:2:nF, [1 2 3]) = [L L + 1 L + nX + 1];


fprintf('mesh generation completed\n\n');


end