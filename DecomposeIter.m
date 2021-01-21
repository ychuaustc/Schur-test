function [verticeAdd, vertice, MC, num] = DecomposeIter(verticeAdd, vertice, MC, num, nV, nv)
%
%   this function compute one decomposition by:
%   1. take one point as the decomposition set (vertice)
%   2. iteratively adding the neighbour (verticeAdd) of current decomposition set (vertice), and add
%   it to vertice to get the new decomposition set
%   3. stop when the size of the current decomposition reaches a threshold: nv
%
%   INPUT:  verticeAdd - vertice to be added in current iteration
%           vertice - the decomposition obtained in current iteration
%           MC - adjacent matrix
%           num - size of vertice before current iteration
%           nV - mesh size
%           nv - expected size of each decomposition
%
%   OUTPUT: verticeAdd - vertice to be added in current iteration
%           vertice - the decomposition obtained in current iteration
%           MC - adjacent matrix
%           num - size of vertice before current iteration


%	compute the neighbour of vertisAdd and update the decomposition size
tempNeighbor = double(double(sum(MC(:, find(verticeAdd)'), 2) > 0) > vertice);
numTemp = num + sum(tempNeighbor);
    

%	do the conditional iteration
if numTemp > num && numTemp <= nv	% if the size has 1. been changed and 2. stays below a given limit, do the iteration
    num = numTemp;
    verticeAdd = tempNeighbor;
    vertice = vertice + verticeAdd;
    [verticeAdd, vertice, MC, num] = DecomposeIter(verticeAdd, vertice, MC, num, nV, nv);
else
end


end