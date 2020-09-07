%%  return a decomposition of given mesh
%%
%%  Input:  n: upper limit of the size of each decomposition (const)
%%          nV: mesh size (const)
%%          B: boundary (column vector)
%%          MC: connection matrix (sparse nV * nV matrix)
%%  Output: DV: list of vertice in the decomposition (tuple of sparse nV * 1 matrix)
%%          DB: list of boundary in the decomposition (tuple of sparse nV * 1 matrix)
%%          DM: list of connection matrix in the decomposition (tuple of sparse nV * nV matrix)
%%          DS: list of subdomain in the decomposition (tuple of sparse nV * 1 matrix)
%%          DE: edges in the decomposition (sparse nV * 1 matrix)
%%          DW: wirebasket in the decomposition (sparse nV * 1 matrix)
%%          numDecompose: number of decompositions (constant)

function [DV, DB, DM, DS, DE, DW, DWB, numDecompose] = Decompose(n, nV, B, MC)

%% initialization
ML = MC; % connection matrix for the left part after one decomposition step
verticeL = ones(nV, 1); % vertice in the left part after one decomposition step
boundL = zeros(nV, 1); % boundary for vertisL
boundL(B) = 1;
column = 1; % step count

%%
fprintf('starting decomposition ...\n');

%%
tic;

%% decomposition step loop
while sum(boundL) > 0
    % pick up a random vertex on Bound, add it to vertice
    boundTemp = find(boundL);
    verticeAdd = zeros(nV, 1);
    vertice = zeros(nV, 1);
    verticeAdd(boundTemp(1)) = 1;
    vertice(boundTemp(1)) = 1;
    num = 1;
    
    % start the decomposition with the vertex added in last step, and save the result to vertice
    [~, vertice, ~, ~] = DecomposeByNeighbor(verticeAdd, vertice, ML, num, n, nV);
    
    % compute the connection matrix for this decomposition (with interface)
    vertisLWOI = verticeL - vertice; % WOL means without interface
    vertisLWOITemp = find(vertisLWOI)';
    MD = ML; % connection matrix for the decomposition
    MD(vertisLWOITemp, :) = 0; %
    MD(:, vertisLWOITemp) = 0; %
    
    % compute the interface between the decomposition and the left part
    MLWOI = ML - MD;
    verticeTemp = find(vertice)';
    InterfaceDTemp = verticeTemp(find(sum(MLWOI(:, verticeTemp), 1))');
    InterfaceD = zeros(nV, 1);
    InterfaceD(InterfaceDTemp) = 1;
    
    % modification of decomposition, connection matrix, interface and boundary
    for i = InterfaceDTemp % if a vertex is in the left part, but only connected with the decomposition (actually the 
                           % interface), then add it to the decomposition
        for j = find(ML(:, i) > InterfaceD)'
            if sum(ML(:, j) > InterfaceD) == 0
                vertice(j) = 1; % exact vertice
                neighborTemp = find(ML(:, j))';
                MD(neighborTemp, j) = 1;
                MD(j, neighborTemp) = 1;  % exact connection matrix for decomposition
            end
        end
    end
    MLWOI = ML - MD;
    verticeTemp = find(vertice)';
    InterfaceDTemp = verticeTemp(find(sum(MLWOI(:, verticeTemp), 1))');
    InterfaceD = zeros(nV, 1);
    InterfaceD(InterfaceDTemp) = 1; % exact interface
    boundD = InterfaceD | (boundL & vertice); % exact boundary for decomposition
    
    % construct decomposition list
    DV{column} = vertice;	% lists of vertice in decomposition
    DB{column} = boundD;	% lists of boundary in decomposition
    DM{column} = MD;     % lists of connection matrix in decomposition
    
    % compute the boundary, connection matrix and vertis for left part
    boundL = (boundL > vertice) | InterfaceD;
    MI = ML;
    vertisWOI = verticeL - InterfaceD;
    vertisWOITemp = find(vertisWOI)';
    MI(vertisWOITemp, :) = 0;
    MI(:, vertisWOITemp) = 0;
    ML = MLWOI + MI; % exact connection matrix for left part
    vertisLWOI = verticeL - vertice;
    verticeL = vertisLWOI | InterfaceD; % exact left part
    
    % update the loop count
    column = column + 1;
end

%% merge small decomposition into large ones
numDecompose = size(DV, 2);
nMerge = ceil(n / 10); % if the size of one decomposition does not reach this limit, we merge it into another one
for i = numDecompose: -1 :1
    if sum(DV{i}) <= nMerge
        for j = 1:numDecompose
            if j ~= i && sum(DV{j} & DV{i}) == 0 && sum(DV{j}) > nMerge
                DV{j} = DV{j} | DV{i};
                DM{j} = DM{j} | DM{i};
                DB{j} = DB{j} | DB{i};
                DV(i) = [];
                DB(i) = [];
                DM(i) = [];
                break;
            end
        end
    end
end
numDecompose = size(DV, 2); % exact # of decomposition

%%
t = toc;

%% check the validation of decomposition
che1 = check1(nV, B, DV, DB, numDecompose);
if che1 == 0
    fprintf('decomposition done\nrunning time for decomposition process: %f s\n\n', t);
else
    fprintf('decomposition failed\n');
    quit(1);
end

%% modification of decomposition list
Bound = zeros(nV, 1);
Bound(B) = 1;
for i = 1:numDecompose
    DS{i} = DV{i} - DB{i}; % exact subdomain list
end

%% find wirebasket points in 2D case
fprintf('computing wirebasket set ...\n');
DSall = zeros(nV, 1);
for i = 1:numDecompose
    DSall = double(DSall | DS{i});
end
dsall = find(DSall);
% wirebasket
DW = ones(nV, 1);
DW(dsall') = 0;
DW(B) = 0;
DW(find(sum(MC(:, dsall'), 2))') = 0; % DW = zeros(nV, 1);
dw = find(DW); % exact wirebasket
% edges
DE = ones(nV, 1);
DE(dsall') = 0;
DE(B) = 0;
DE(dw') = 0;
for i = 1:numDecompose
    DB{i} = double(DB{i} > Bound);
    DB{i} = double(DB{i} > DW); % exact boundary list
end
%
DWB = zeros(nV, 1);
DWB(find(sum(MC(:, dw'), 2))') = 1;
DWB = double(DE & DWB); % boundary of wirebasket

%% check the validation of wirebasket
che2 = check2(nV, DS, B, DE, DW, numDecompose);
if che2 == 0
    fprintf('wirebasket set computed\n\n');
else
    fprintf('wirebasket computation failed\n');
    quit(1);
end

%% decomposition information
IntersectionWOB = zeros(nV, 1); % union of boundaries of all decompositions
for i = 1:numDecompose
	IntersectionWOB = IntersectionWOB | DB{i};
end
numIntersectionWOB = sum(IntersectionWOB);
for i = 1:numDecompose
    numDS(1, i) = sum(DS{i});
end
minSubdomain = min(numDS);
maxSubdomain = max(numDS);
meanSubdomain = floor(sum(numDS) / numDecompose);

fprintf('decomposition information:\n');
fprintf('number of decomposition: %d \n', numDecompose);
fprintf('minimum size of subdomain in decomposition: %d \n', minSubdomain);
fprintf('maximum size of subdomain in decomposition: %d \n', maxSubdomain);
fprintf('mean size of subdomain in decomposition: %d \n', meanSubdomain);
fprintf('size of intersection: %d \n\n\n\n', numIntersectionWOB);

end