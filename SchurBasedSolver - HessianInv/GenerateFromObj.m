function [Vertex, Face] = GenerateFromObj(Filename)
%
%   this function read obj file and generate the triangulation by calling the readObj function
%
%   INPUT:  Filename - the obj file
%
%   OUTPUT: Vertex - mesh vertices
%           Face - mesh faces


fprintf('mesh generating...\n');


[Vertex, Face] = readObj(Filename);


fprintf('mesh generation completed\n\n');


end