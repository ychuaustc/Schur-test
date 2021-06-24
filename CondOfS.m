function [condS1, condS2] = CondOfS(MSS, MSE, MWW, MWE, MEE, MEE_W, numDecompose)
%
%   compute condition number of Schur matrix
%
%   INPUT:  
%
%   OUTPUT: 


S = MEE;
for k = 1:size(MEE, 2)
    for i = 1:numDecompose
        V{i} = MSS{i} \ MSE{i}(:, k);
    end
    for i = 1:numDecompose
        S(:, k) = S(:, k) - MSE{i}' * V{i};
    end
end
S = S - MWE' * (MWW \ MWE);
condS1 = condest(S);

for k = 1:size(S, 2)
    S(:, k) = Preconditioner1(MWW, MWE, MEE_W, S(:, k));
end
condS2 = condest(S);


end