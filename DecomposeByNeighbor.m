%%  return the result (vertice, connection matrix, etc.) for one decomposition step
%%
%%  for one decomposition step, the code processes as follows:
%%      pick up one vertex on the boundary of the vertice set, add it to the decomposition, when the size of the 
%%      decomposition is less than the given limit, do the iteration:
%%          find the neighbor of all added vertis in last iteration, and add them to the decomposition

%%  Input:
%%      verticeAdd: vertice added in last iteration (nV * 1 sparse matrix)
%%      vertice: vertice gained after last iteration (nV * 1 sparse matrix)
%%      M: connection matrix for vertice (nV * nV sparse matrix; stay unchanged in all iterations, used only for 
%%      the computation of the neighborhood)
%%      num: size of vertice after last iteration (constant)
%%      n: upper limit of the size of one decomposition (constant)
%%      nV: mesh size
%%  Output:
%%      verticeAdd: vertice added in current iteration (nV * 1 sparse matrix)
%%      vertice: vertice gained after current iteration (nV * 1 sparse matrix)
%%      M: connection matrix for vertice (nV * nV sparse matrix; stay unchanged in all iterations, used only for computing the 
%%         neighborhood)
%%      num: size of vertice after current iteration (constant)

function [verticeAdd, vertice, M, num] = DecomposeByNeighbor(verticeAdd, vertice, M, num, n, nV)

% compute neighbor for vertisAdd and update the decomposition size
tempNeighbor = double(double(sum(M(:, find(verticeAdd)'), 2) > 0) > vertice);
numTemp = num + sum(tempNeighbor);
    
% do the conditional iteration
if numTemp > num && numTemp <= n % if the size has been changed and still below the given limit, do the iteration
    num = numTemp;
    verticeAdd = tempNeighbor;
    vertice = vertice + verticeAdd;
    [verticeAdd, vertice, M, num] = DecomposeByNeighbor(verticeAdd, vertice, M, num, n, nV);
else
end

end