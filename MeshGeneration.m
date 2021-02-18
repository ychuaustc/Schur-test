function [Vertex, Face] = MeshGeneration(meshType, nV, fileName)
%
%
%   this function generates triangulations upon mesh type and mesh size
%
%   INPUT:  meshType - the type of the triangulation and the manner to obtain it
%           nV - mesh size
%           fileName - name of obj file
%
%   OUTPUT: Vertex - mesh vertices
%           Face - mesh faces


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