function [Vertex, Face] = GenerateFromObj(Filename)
%
%   this function reads obj file and generate the triangulation by calling the readObj function
%
%   INPUT:  Filename - the obj file name
%
%   OUTPUT: Vertex - mesh vertices
%           Face - mesh faces


[Vertex, Face] = readObj(Filename);


end