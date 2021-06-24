function [fwe] = FindFforWE(FinWE, dweInd)

fwe = FinWE;
nwe = size(fwe, 1);
for i = 1:3
    for j = 1:nwe
        fwe(j, i) = find(dweInd == fwe(j, i));
    end
end