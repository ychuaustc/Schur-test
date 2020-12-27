%%  check the validation of the decomposition
%%
%%  Input:
%%          nV:             mesh size
%%          B:              boundary
%%          DV:             list of vertice in the decomposition
%%          DB:             list of boundary in the decomposition
%%          numDecompose:	number of decomposition
%%          t:              running time
%%  Output: 

function [] = DecomposeUTCheck1(nV, B, DV, DB, numDecompose, t)

%%
test1 = zeros(nV, 1); %	union of intersections of all decompositions and boundary of the mesh
test2 = zeros(nV, 1); %	union of boundaries of all decompositions
testB = zeros(nV, 1);
testB(B) = 1;
for i = 1:numDecompose - 1
    for j = i+1:numDecompose
        test1 = test1 | (DV{i} & DV{j});
    end
end
test1 = test1 | testB;
for i = 1:numDecompose
	test2 = test2 | DB{i};
end
che1 = sum(test1 - test2);

%%
if che1 == 0
    fprintf('decomposition done\nrunning time for decomposition process: %f s\n\n', t);
else
    fprintf('decomposition failed\n');
    quit(1);
end

%%
end