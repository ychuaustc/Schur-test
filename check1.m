%%  check the validation of the decomposition
%%
%%  Input:  nV: mesh size (const)
%%          B: boundary (column vector)
%%          DV: list of vertice in the decomposition (tuple of sparse nV * 1 matrix)
%%          DB: list of boundary in the decomposition (tuple of sparse nV * 1 matrix)
%%          numDecompose: number of decomposition (constant)
%%  Output: che1: check result (constant)

function [che1] = check1(nV, B, DV, DB, numDecompose)

%%
test1 = zeros(nV, 1); % union of intersections of all decompositions and boundary of the mesh
test2 = zeros(nV, 1); % union of boundaries of all decompositions
testB = zeros(nV, 1);
testB(B) = 1;
for i = 1:numDecompose-1
    for j = i+1:numDecompose
        test1 = test1 | (DV{i} & DV{j});
    end
end
test1 = test1 | testB;
for i = 1:numDecompose
	test2 = test2 | DB{i};
end
che1 = sum(test1 - test2);

end