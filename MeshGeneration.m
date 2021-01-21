function [Vertex, Face] = MeshGeneration(meshType, nV)
%
%
%   this function generates triangulations upon the mesh type and the mesh size
%
%   INPUT:  Vertex - mesh vertices
%           Face - mesh faces
%
%   OUTPUT: meshType - the type of the triangulation and the manner to get it
%           nV - mesh size


switch meshType
    case 1
        [Vertex, Face] = GenerateST(nV);
    case 2
        [Vertex, Face] = GenerateUT(nV);
    case 3
        [Vertex, Face] = GenerateFromObj('NewtonParam\ball');
    otherwise
        quit(1);
end


end