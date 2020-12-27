%%  return the result (vertice, connection matrix, etc.) for one decomposition step
%%
%%  for one decomposition step:
%%      pick up one vertex on the boundary of the vertice set, add it to the decomposition, until the size of the 
%%      decomposition is less than a given limit, do the iteration:
%%          find the neighbor of all added vertice in last iteration, and add them to the decomposition
%%
%%  Input:
%%      verticeAdd: vertice added in last iteration
%%      vertice:    vertice gained after last iteration
%%      M:          connection matrix for vertice, stay unchanged during all iterations
%%      num:        size of vertice after last iteration
%%      nV:         mesh size
%%      nv:         upper limit of the size of one decomposition
%%  Output:
%%      verticeAdd:	vertice added in current iteration
%%      vertice:    vertice gained after current iteration
%%      M:          connection matrix for vertice, stay unchanged during all iterations
%%      num:        size of vertice after current iteration

function [verticeAdd, vertice, M, num] = DecomposeByNeighbor(verticeAdd, vertice, M, num, nV, nv)

%	compute neighbor for vertisAdd and update the decomposition size
tempNeighbor = double(double(sum(M(:, find(verticeAdd)'), 2) > 0) > vertice);
numTemp = num + sum(tempNeighbor);
    
%	do the conditional iteration
if numTemp > num && numTemp <= nv %	if the size has been changed but still stays below a given limit, do the iteration
    num = numTemp;
    verticeAdd = tempNeighbor;
    vertice = vertice + verticeAdd;
    [verticeAdd, vertice, M, num] = DecomposeByNeighbor(verticeAdd, vertice, M, num, nV, nv);
else
end

%%
end