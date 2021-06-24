function [DS, DE, DW, numDecompose] = Decomp1()

numDecompose = 1;

DS{1} = [0 0 0 0 0 0 1 1]';
DW = [1 1 0 0 0 0 0 0 0]';
DE = [0 0 1 1 1 1 0 0]';

end