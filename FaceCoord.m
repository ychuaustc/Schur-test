%%  return coordinate of triangles in Face
%%
%%  Input:
%%          Vertex:     vertice matrix
%%          Face:       face matrix
%%  Output:
%%          FX:         x coordinate for each triangle in Face
%%          FY:         y coordinate for each triangle in Face

function [FX, FY] = FaceCoord(Vertex, Face)

%%
%
VX = Vertex(:, 1);
VY = Vertex(:, 2);
FX = VX(Face);
FY = VY(Face);

%%
end