%%  
function [] = testfunc(Vertex, Face, DS, DE, DW, DWB, nV, nF, numDecompose)
%%
M1 = findM(Vertex, Face);
[MSS, MSE, MWW, MWE, MEE, dsInd, dwInd, deInd] = SchurSystemM(M1, DS, DE, DW, numDecompose);
%%
FNotInWWB = zeros(nF, 1);
Face1 = Face(:, 1);
Face2 = Face(:, 2);
Face3 = Face(:, 3);
dwwbInd = find(DW | DWB)';
notindwwbInd = setdiff(1:nV, dwwbInd);
for i = notindwwbInd
    fnotinwwb1 = find(Face1 == i)';
    FNotInWWB(fnotinwwb1) = 1;
    fnotinwwb2 = find(Face2 == i)';
    FNotInWWB(fnotinwwb2) = 1;
    fnotinwwb3 = find(Face3 == i)';
    FNotInWWB(fnotinwwb3) = 1;
end
fnotinwInd = find(FNotInWWB)';
finwInd = setdiff(1:nF, fnotinwInd);
Vertex1 = Vertex(dwwbInd, :);
Face1 = Face(finwInd, :);
MEE_W = findM(Vertex, Face1);

end