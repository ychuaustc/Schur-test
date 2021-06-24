function [FinWE, finweInd] = FindFinWE(Face, dweInd, nV)

notinwe = setdiff(1:nV, dweInd);

F1 = Face(:, 1);
F1(ismember(F1, notinwe))=0;
F2 = Face(:, 2);
F2(ismember(F2, notinwe))=0;
F3 = Face(:, 3);
F3(ismember(F3, notinwe))=0;

finweInd = find(F1 & F2 & F3)';
FinWE = Face(finweInd, :);

end