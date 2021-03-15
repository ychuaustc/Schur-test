function [C] = ArapSystemC(MThetaLV)
%
%   this function compute the right hand side of the ARAP system
%
%   INPUT:  MThetaLV - assign the value of cot(theta) * L(V) to the adjacent matrix
%
%   OUTPUT: C - the right hand side of the ARAP system


C = sum(MThetaLV, 2) - sum(MThetaLV, 1)';


end