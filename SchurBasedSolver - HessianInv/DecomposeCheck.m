function [] = DecomposeCheck(nV, B, DV, DB, numDecompose, t)
%
%   this function checks the validation of the decomosition set computed: does the union of the
%   intersections of all decomposition sets match the union of the boundaries of all decomposition
%   sets
%
%   INPUT:  nV - mesh size
%           B - boundary
%           DV - decomposition list
%           DB - subdomain boundary list
%           numDecompose - decomposition number
%           t - running time
%
%   OUTPUT: 


test1 = zeros(nV, 1);
test2 = zeros(nV, 1);
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
if che1 == 0
    fprintf('decomposition done\nrunning time for the decomposition process: %f s\n\n', t);
else
    fprintf('decomposition failed\n');
    quit(1);
end


end