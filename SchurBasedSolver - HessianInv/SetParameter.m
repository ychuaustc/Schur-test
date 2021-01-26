function [meshType, nV, numDecomposeTemp, nv, wirebasketType, solverType, epsArap, epsSchur, fileName] = SetParameter()
%
%   this function sets up the global parameters, such as mesh size, decomposition number, wirebasket
%   type, convergence error and so on
%
%   INPUT:  
%
%   OUTPUT: nV - mesh size
%           numDecomposeTemp - expected decomposition number
%           nv - expected size of each decomposition
%           meshType - the type of the triangulation and the manner to get it
%           wirebasketType - the type of the wirebasket set
%           solverType - the type of the solver
%           epsArap - the convergence error of the local-global method
%           epsSchur - the convergence error of the conjugate method
                   

meshType = listdlg('PromptString', 'Trianulation generating type', ...
                   'SelectionMode','single', ...
                   'liststring', {'Generate structured triangulation directly', ...
                                  'Generate unstructured delaunay triangulation', ...
                                  'Read triangulation from file'});


fileName = 0;
switch meshType
    case {1, 2}
        m1 = listdlg('PromptString', 'Mesh size', 'SelectionMode','single', ...
                     'liststring', {'289', '1089', '4225'}) + 3;
        nV = (2^m1 + 1)^2;
    case 3
        fileName = listdlg('PromptString', 'Obj file name', ...
                           'SelectionMode','single', ...
                           'liststring', {'ball', 'camelhead_slim', 'face'});
        switch fileName
            case 1
                [V, ~, ~, ~] = readObj('NewtonParam/ball');
            case 2
                [V, ~, ~, ~] = readObj('NewtonParam/camelhead_slim');
            case 3
                [V, ~, ~, ~] = readObj('NewtonParam/face');
            otherwise
                quit(1);
        end
        nV = size(V, 1);
    otherwise
        quit(1);
end


if nV <= 4000
	m2 = listdlg('PromptString', 'Decomposition number', 'SelectionMode','single', 'liststring', {'4'});
else
    m2 = listdlg('PromptString', 'Decomposition number', 'SelectionMode','single', 'liststring', {'4', '16'});
end
numDecomposeTemp = 2^(2 * m2);


energyType = listdlg('PromptString', 'Energy type', 'SelectionMode','single', ...
                     'liststring', {'ARAP', 'Symmetric Dirichlet', 'SARAP'});
switch energyType
    case 1
        addpath(genpath('arap'));
    case 2
        addpath(genpath('symmd'));
    case 3
        addpath(genpath('sarap'));
    otherwise
        quit(1);
end


switch meshType
    case {1, 2}
        nv = (2^(m1 - m2) + 1)^2;
    case 3
        nv = floor(nV / numDecomposeTemp);
    otherwise
        quit(1);
end


wirebasketType = 1;
solverType = listdlg('PromptString', 'Solver type', 'SelectionMode','single', ...
                     'liststring', {'Direct solver', 'Conjugate solver', 'Preconditioned conjugate solver'});


eps = floor(1.5 * log10(1.0 / nV));
epsArap = 10^eps;
epsSchur = 10^eps;


end