%%  return information of linear system
%%
%%  Input:  Vertex: vertice matrix (sparse matrix, vertis size * 3)
%%          B / I: boundary/ inner part of mesh (column vector, which denotes the position in vertice)
%%          nB / nI: size of B / I
%%          MC: connection matrix (sparse nV * nV matrix)
%%          DS: list of subdomain of the decomposition (tuple of sparse nV * 1 matrix)
%%          DE: edges in the decomposition (sparse nV * 1 matrix)
%%          DW: wirebasket of the decomposition (sparse nV * 1 matrix)
%%          numDecompose: number of decomposition (constant)
%%  Output: //

function [M, Cx, Cy, dsInd, dwInd, deInd, debInd, dewbInd, MS, MSE, MW, MWE, ME, MB, MSB, MWB, MWWB, CSx, CSy, CWx, CWy, CEx, CEy] ...
         = Lsystem(Vertex, B, nB, I, nI, MC, DS, DE, DW, DB, DWB, numDecompose)

%%
fprintf('setting linear system ...\n');

%%
vertex = Vertex;

%% compute the connection matrix for inner vertex
MIn = MC(I, :); % the modified connection matrix, each row of the matrix represents the connection information of an inner vertex

%% fix the boundary on unit square [0, 1] * [0, 1] * 0
vertex(:, 3) = 0;
vertex(B, 1:2) = SquareBoundary(nB);

%% compute coefficient matrix
M = MIn;
E = sparse(1:nI, 1:nI, sum(M, 2), nI, nI);
C = zeros(nI, 1);

Cx = C - M(:, B) * vertex(B, 1);
Cy = C - M(:, B) * vertex(B, 2);

M(:, B) = [];
M = M - E;
    
%% reform the coefficient matrix upon decomposition

% compute index for each matrix block
% subdomain
for i = 1:numDecompose
    DSInd{i} = DS{i}(I);
    dsInd{i} = find(DSInd{i}); % subdomain index in inner part
end
% wirebasket
DWInd = DW(I);
dwInd = find(DWInd); % wirebasket id in inner part
% edges
DEInd = DE(I);
deInd = find(DEInd); % edge index in inner part
for i = 1:numDecompose
    DBInd{i} = DB{i}(I);
    dbInd{i} = find(DBInd{i}); % boundary index in inner part
end
% debInd
for i = 1:numDecompose
    debInd{i} = find(DBInd{i}(deInd));
end
% debwInd
DWBInd = DWB(I);
dwbInd = find(DWBInd);
dewbInd = find(DWBInd(deInd));

% define matrix blocks with index computed
for i = 1:numDecompose
    MS{i} = M(dsInd{i}, dsInd{i}); % matrix block for subdomain
end
for i = 1:numDecompose
    MSE{i} = M(dsInd{i}, deInd); % matrix block for subdomain intersect edges
end
MW = M(dwInd, dwInd);
MWE = M(dwInd, deInd); % matrix block for wirebasket intersect edges
ME = M(deInd, deInd);
for i = 1:numDecompose
    MB{i} = M(dbInd{i}, dbInd{i}); % matrix block for boundary
end
for i = 1:numDecompose
    MSB{i} = M(dsInd{i}, dbInd{i}); % matrix block for boundary
end
MWB = M(dwbInd, dwbInd);
MWWB = M(dwInd, dwbInd);
for i = 1:numDecompose
    CSx{i} = Cx(dsInd{i});
    CSy{i} = Cy(dsInd{i});
end
CWx = Cx(dwInd);
CWy = Cy(dwInd);
CEx = Cx(deInd);
CEy = Cy(deInd);

%%
fprintf('linear system set\n\n\n\n');

end