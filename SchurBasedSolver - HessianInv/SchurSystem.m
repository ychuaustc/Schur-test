%%	construct sparse system MX = C with given decomposition
%%
%%	Input:
%%          M:              M matrix
%%          CT:             C vector (tuple)
%%          DS:             list of subdomain in the decomposition
%%          DB:             list of boundary in the decomposition
%%          DE:             edges in the decomposition
%%          DW:             wirebasket in the decomposition
%%          DWB:            boundary of wirebasket in the decomposition           
%%          numDecompose:   number of decompositions
%%	Output:
%%          MSS:            matrix block for subdomain * subdomain (tuple)
%%          MSE:            matrix block for subdomain * edge (tuple)
%%          MWW:            matrix block for wirebasket * wirebasket
%%          MWE:            matrix block for wirebasket * edge
%%          MEE:            matrix block for edge * edge
%%          CST:            matrix block for subdomain (2d tuple)
%%          CWT:            matrix block for wirebasket (tuple)
%%          S:              Schur matrix
%%          bT:             Schur vector (tuple)
%%          dsInd:          index for subdomain (tuple)
%%          dwInd:          index for wirebasket (tuple)
%%          deInd:          index for edge (tuple)

function [MSS, MSE, MWW, MWE, MEE, CST, CWT, S, bT, dsInd, dwInd, deInd] = SchurSystem(M, CT, DS, DB, DE, DW, DWB, numDecompose)

%%
fprintf('constructing Schur system ...\n');

%%
tic;

%% compute decomposition index
%  subdomain
for i = 1:numDecompose
    dsInd{i} = find(DS{i}); % subdomain index
end
%  wirebasket
dwInd = find(DW); % wirebasket index
%  edge
deInd = find(DE); % edge index
%  boundary
for i = 1:numDecompose
    dbInd{i} = find(DB{i}); % boundary index
end
dwbInd = find(DWB); % wirebasket boundary index

%%  define matrix blocks for M
for i = 1:numDecompose
    MSS{i} = M(dsInd{i}, dsInd{i}); % matrix block for subdomain * subdomain
end
for i = 1:numDecompose
    MSE{i} = M(dsInd{i}, deInd); % matrix block for subdomain * edge
end
MWW = M(dwInd, dwInd); % matrix block for wirebasket * wirebasket
MWE = M(dwInd, deInd); % matrix block for wirebasket * edge
MEE = M(deInd, deInd); % matrix block for edge * edge

%%  define matrix blocks for CT
nC = size(CT, 2); % size of tuple CT
for i = 1:nC
    for j = 1:numDecompose
        CSTemp{j} = CT{i}(dsInd{j});
    end
    CST{i} = CSTemp;
    CWT{i} = CT{i}(dwInd);
    CET{i} = CT{i}(deInd);
end

%%  compute Schur matrix S and Schur vector b
%   Schur matrix S
S = MEE;
for i = 1:numDecompose
    S = S - MSE{i}' * (MSS{i} \ MSE{i});
end
S = S - MWE' * (MWW \ MWE);
%	Schur vector b
for i = 1:nC
    bT{i} = CET{i};
    for j = 1:numDecompose
        bT{i} = bT{i} - MSE{j}' * (MSS{j} \ CST{i}{j});
    end
    bT{i} = bT{i} - MWE' * (MWW \ CWT{i});
end

%%
t = toc;

%%
fprintf('Schur system constructed\nrunning time for Schur system construction: %f s\n\n\n\n', t);

%%
end