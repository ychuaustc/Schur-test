function [X, t] = DirectSolver(M, C)
%
%   this function solve MX = C directly
%
%   INPUT:  M - the coefficient matrix
%           C - the right hand side of the system
%
%   OUTPUT: X - the solution
%           t - running time


tic;

X = M \ C;

t = toc;


end