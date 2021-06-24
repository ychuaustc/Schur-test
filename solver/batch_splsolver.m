classdef batch_splsolver < handle
    properties(SetAccess = protected)
        handle = -1;
        A = {};  % matrix, nonzero pattern
    end
    
    methods(Static)
        function i = solverid(solvername)
            switch lower(solvername)
                case {'llt'}
                    i = 0;
                case {'ldlt'}
                    i = 1;
                case {'lu'}
                    i = 2;
                otherwise
                    error('unknown solver %s\n', solvername);
            end
        end
        
        function x = fullsolve(A, b, solvername)
            if nargin<3, solvername = 'lu'; end
%             tic
            [x, f] = batch_splsolver_imp(5, A, b, batch_splsolver.solverid(solvername));
%             fprintf('splsolver full      = %fs\n', toc);
            assert(all(f), 'linear solver failed\n');
        end
    end

    methods
        function s = batch_splsolver(A, solvername)   %A is a cell & pattern only
            if nargin<2, solvername = 'lu'; end
            s.A = A;
%             tic
            s.handle = batch_splsolver_imp(1, A, batch_splsolver.solverid(solvername));
%             fprintf('splsolver    symf   = %fs\n', toc);
        end

        function refactorize(s, nonzeros)
            % TODO: checke nonzeros is a vector
%             tic
            f = batch_splsolver_imp(2, s.handle, s.A, nonzeros);
%             fprintf('splsolver    numf   = %fs\n', toc);
            assert(all(f), 'linear solver failed with factorization');
        end
        
        function x = solve(s, b)
%             tic
            x = batch_splsolver_imp(3, s.handle, b);
%             fprintf('splsolver    solve  = %fs\n', toc);
        end
        
        function x = mldivide(s, b)
%             tic
            x = batch_splsolver_imp(3, s.handle, b);
%             fprintf('splsolver    solve  = %fs\n', toc);
        end
        
        function [x, f] = refactor_solve(s, nonzeros, b)
%             tic
            [x, f] = batch_splsolver_imp(4, s.handle, s.A, nonzeros, b);
%             fprintf('splsolver    f&s    = %fs\n', toc);
            assert(all(f), 'linear solver failed');
        end
        
        function delete(s)
            if s.handle ~= -1
%                 tic
                batch_splsolver_imp(6, s.handle);
%                 fprintf('splsolver    clean up   = %fs\n', toc);
            end
        end
    end
end