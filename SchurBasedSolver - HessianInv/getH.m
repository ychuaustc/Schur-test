function [MH] = getH(meshType, nV, fileName, numDecompose)


switch meshType
    case 1
        switch nV
            case 289
                M = load('H/s_289_4');
            case 1089
                M = load('H/s_1089_4');
            case 4225
                if numDecompose < 10
                    M = load('H/s_4225_4');
                else
                    M = load('H/s_4225_16');
                end
            otherwise
                quit(1);
        end
    case 2
        switch nV
            case 289
                M = load('H/u_289_4');
            case 1089
                M = load('H/u_1089_4');
            case 4225
                if numDecompose < 10
                    M = load('H/u_4225_4');
                else
                    M = load('H/u_4225_16');
                end
            otherwise
                quit(1);
        end
    case 3
        switch fileName
            case 1
                M = load('H/ball_4');
            case 2
                if numDecompose < 10
                    M = load('H/camel_4');
                else
                    M = load('H/camel_16');
                end
            case 3
                M = load('H/face_4');
            otherwise
                quit(1);
        end
    otherwise
        quit(1);
end


MH = M.H{3};


end