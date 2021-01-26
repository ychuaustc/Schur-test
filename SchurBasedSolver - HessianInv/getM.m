function [MH] = getM(meshType, nV, fileName)


switch meshType
    case 1
        switch nV
            case 289
                M = load('M/s_289');
            case 1089
                M = load('M/s_1089');
            case 4225
                M = load('M/s_4225');
            otherwise
                quit(1);
        end
    case 2
        switch nV
            case 289
                M = load('M/u_289');
            case 1089
                M = load('M/u_1089');
            case 4225
                M = load('M/u_4225');
            otherwise
                quit(1);
        end
    case 3
        switch fileName
            case 1
                M = load('M/ball');
            case 2
                M = load('M/camel');
            case 3
                M = load('M/face');
            otherwise
                quit(1);
        end
    otherwise
        quit(1);
end


MH = M.MH;


end