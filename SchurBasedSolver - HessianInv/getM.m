function [MH] = getM(meshType, nV, fileName)


switch meshType
    case 1
        switch nV
            case 289
                M = load('Hessian/MS_289');
            case 1089
                M = load('Hessian/MS_1089');
            case 4225
                M = load('Hessian/MS_4225');
            otherwise
                quit(1);
        end
    case 2
        switch nV
            case 289
                M = load('Hessian/MU_289');
            case 1089
                M = load('Hessian/MU_1089');
            case 4225
                M = load('Hessian/MU_4225');
            otherwise
                quit(1);
        end
    case 3
        switch fileName
            case 1
                M = load('Hessian/M_Ball');
            case 2
                M = load('Hessian/M_Camel');
            case 3
                M = load('Hessian/M_Face');
            otherwise
                quit(1);
        end
    otherwise
        quit(1);
end


MH = M.MH;


end