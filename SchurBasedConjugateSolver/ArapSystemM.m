function [M] = ArapSystemM(MCotTheta, nV)
%
%   this function computes the coefficient matrix for the ARAP system in the global step
%
%   INPUT:  MCotTheta - the cotangent value of the adjacent matrix
%           nV - mesh size
%
%   OUTPUT: M - the coefficient matrix


M = sparse(1:nV, 1:nV, sum(MCotTheta + MCotTheta', 2), nV, nV) - (MCotTheta + MCotTheta');


end