%%  compute the L1-error of two vertice set
%%
%%  Input:  vec1:	vector 1 of size n
%%          vec2:	vector 2 of size n
%%          n:      vector size
%%  Output: err:    L1-error

function [err] = L1error(vec1, vec2, n)

err = sum(abs(vec1 - vec2)) / n;

end