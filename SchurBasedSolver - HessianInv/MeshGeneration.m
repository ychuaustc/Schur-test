function [Vertex, Face] = MeshGeneration(meshType, nV, fileName)
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
        switch fileName
            case 1
                [Vertex, Face] = GenerateFromObj('NewtonParam/ball');
            case 2
                [Vertex, Face] = GenerateFromObj('NewtonParam/camelhead_slim');
            case 3
                [Vertex, Face] = GenerateFromObj('NewtonParam/face');
            otherwise
                quit(1);
        end
    otherwise
        quit(1);
end


end