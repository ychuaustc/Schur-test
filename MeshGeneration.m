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
                [Vertex, Face] = GenerateFromObj('NewtonParam/bunny');
            case 3
                [Vertex, Face] = GenerateFromObj('NewtonParam/camel');
            case 4
                [Vertex, Face] = GenerateFromObj('NewtonParam/camelhead45k');
            case 5
                [Vertex, Face] = GenerateFromObj('NewtonParam/camelhead181k');
            case 6
                [Vertex, Face] = GenerateFromObj('NewtonParam/botijo_parametrized_nonlin_cut');
            case 7
                [Vertex, Face] = GenerateFromObj('NewtonParam/bozbezbozzel100K_parametrized_nonlin_cut');
            case 8
                [Vertex, Face] = GenerateFromObj('NewtonParam/blade_parametrized_nonlin_cut');
            case 9
                [Vertex, Face] = GenerateFromObj('NewtonParam/amphora_parametrized_nonlin_cut');
            case 10
                [Vertex, Face] = GenerateFromObj('NewtonParam/aircraft_cut');
            case 11
                [Vertex, Face] = GenerateFromObj('NewtonParam/face');
            otherwise
                quit(1);
        end
    otherwise
        quit(1);
end


end