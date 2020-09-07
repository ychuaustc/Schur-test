%%  return the wirebasket set of a decomposition
%%
%%  Input:
%%      verticeAdd: 
%%      vertice: 
%%      M: 
%%      B: 
%%      nV: mesh size
%%  Output:
%%      verticeAdd: 
%%      vertice: 
%%      M:

function [verticeAdd, vertice, M] = WirebasketEbyNeighbor(verticeAdd, vertice, M, nV)

tempNeighbor = double(double(sum(M(:, find(verticeAdd)'), 2) > 0) > vertice);
    
% do the conditional iteration
if  sum(tempNeighbor) ~= 0 % 
    verticeAdd = tempNeighbor;
    vertice = vertice + verticeAdd;
    [verticeAdd, vertice, M] = WirebasketEbyNeighbor(verticeAdd, vertice, M, nV);
else
end

end