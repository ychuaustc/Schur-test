function [FX, FY] = FaceCoord(Vertex, Face)
%
%   this function assigns coordinates for the face matrix 
%
%   INPUT:  Vertex - mesh vertices
%           Face - mesh faces
%
%   OUTPUT: FX - x-coordinates for the face matrix
%           FY - y-coordinates for the face matrix


VX = Vertex(:, 1);
VY = Vertex(:, 2);
FX = VX(Face);
FY = VY(Face);


end