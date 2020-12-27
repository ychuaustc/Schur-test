%%  return a decomposition of given mesh
%%
%%  Input:
%%          nV:             mesh size
%%          nv:             upper limit of the size of each decomposition
%%          B:              boundary
%%          MC:             connection matrix
%%  Output:
%%          DV:
%%          DS:             list of subdomain in the decomposition
%%          DB:             list of boundary in the decomposition
%%          DE:             edges in the decomposition
%%          DW:             wirebasket in the decomposition
%%          DWB:            boundary of wirebasket in the decomposition
%%          numDS:          size of each decomposition
%%          numDecompose:   number of decompositions

function [DV, DS, DB, DE, DW, DWB, numDS, numDecompose] = Decompose(nV, nv, B, MC)

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

%%	decomposition step loop
while sum(boundL) > 0
    %	pick up a random vertex on Bound, add it to vertice
    boundTemp = find(boundL);
    verticeAdd = zeros(nV, 1);
    vertice = zeros(nV, 1);
    verticeAdd(boundTemp(1)) = 1;
    vertice(boundTemp(1)) = 1;
    num = 1;
    %   do the decomposition with the vertex added in last step, and save the result to vertice
    [~, vertice, ~, ~] = DecomposeByNeighbor(verticeAdd, vertice, ML, num, nV, nv);
    %	compute the connection matrix for this decomposition (with interface)
    verticeLWOI = verticeL - vertice; %	WOL means without interface
    verticeLWOITemp = find(verticeLWOI)';
    MD = ML; %	connection matrix for the decomposition
    MD(verticeLWOITemp, :) = 0;
    MD(:, verticeLWOITemp) = 0;
    %   compute the interface between the decomposition and the left part
    MLWOI = ML - MD;
    verticeTemp = find(vertice)';
    InterfaceDTemp = verticeTemp(find(sum(MLWOI(:, verticeTemp), 1))');
    %	modification of decomposition, connection matrix, interface and boundary: if a vertex is in the left part, but only 
    %	connected with the decomposition (the interface), we add it to the decomposition
    %
    interfaceDNeighborTemp = find(sum(ML(:, InterfaceDTemp), 2))';
    interfaceDNeighborComp = ones(nV, 1);
    interfaceDNeighborComp(InterfaceDTemp) = 0;
    interfaceDNeighborCompTemp = find(interfaceDNeighborComp)';
    NeighborToAddTemp = interfaceDNeighborTemp(find(sum(ML(interfaceDNeighborCompTemp, interfaceDNeighborTemp), 1) == 0)');
    %
    NeighborToAdd = zeros(nV, 1);
    NeighborToAdd(NeighborToAddTemp) = 1;
    NeighborToAddTemp = find(double(NeighborToAdd > vertice))';
    %
    vertice(NeighborToAddTemp) = 1; %	exact decomposition
    %
%     MD(:, NeighborToAddTemp) = ML(:, NeighborToAddTemp);
%     MD(NeighborToAddTemp, :) = ML(NeighborToAddTemp, :); % exact connection matrix
    verticeLWOI = verticeL - vertice;
    verticeLWOITemp = find(verticeLWOI)';
    MD = ML;
    MD(verticeLWOITemp, :) = 0;
    MD(:, verticeLWOITemp) = 0; %	exact connection matrix
    %
    MLWOI = ML - MD;
    verticeTemp = find(vertice)';
    InterfaceDTemp = verticeTemp(find(sum(MLWOI(:, verticeTemp), 1))');
    InterfaceD = zeros(nV, 1);
    InterfaceD(InterfaceDTemp) = 1; %	exact interface
    boundD = InterfaceD | (boundL & vertice); %	exact boundary for decomposition
    % construct decomposition list
    DV{column} = vertice; %	lists of vertice in decomposition
    DB{column} = boundD;  %	lists of boundary in decomposition
    DM{column} = MD;      %	lists of connection matrix in decomposition
    %	compute the boundary, connection matrix and vertice for left part
    boundL = (boundL > vertice) | InterfaceD; %	exact boundary for left part
    MI = ML;
    verticeWOI = verticeL - InterfaceD;
    verticeWOITemp = find(verticeWOI)';
    MI(verticeWOITemp, :) = 0;
    MI(:, verticeWOITemp) = 0;
    ML = MLWOI + MI; %	exact connection matrix for left part
    verticeLWOI = verticeL - vertice;
    verticeL = verticeLWOI | InterfaceD; %	exact left part
    %	update the loop count
    column = column + 1;
end

%%	merge small decomposition into large ones
numDecompose = size(DV, 2);
nMerge = ceil(nv / 5); %	if the size of one decomposition does not reach this limit, we merge it into another one
numDecomposeTemp = numDecompose;
for i = numDecompose: -1 :1
    if sum(DV{i}) <= nMerge
        for j = 1:numDecomposeTemp
            if j ~= i && sum(DV{j} & DV{i}) == 0 && sum(DV{j}) > nMerge && sum(DV{j}) < 4 * nMerge
                DV{j} = DV{j} | DV{i};
                DM{j} = DM{j} | DM{i};
                DB{j} = DB{j} | DB{i};
                DV(i) = [];
                DB(i) = [];
                DM(i) = [];
                numDecomposeTemp = numDecomposeTemp - 1;
                break;
            end
        end
    end
end
numDecompose = size(DV, 2); %	exact # of decomposition

%%
t = toc;

%%	check the validation of decomposition
DecomposeUTCheck1(nV, B, DV, DB, numDecompose, t);

%%	modification of decomposition list
for i = 1:numDecompose
    DS{i} = DV{i} - DB{i}; %	exact list of subdomain
end

%%	find wirebasket set
fprintf('computing wirebasket set ...\n');
DSall = zeros(nV, 1);
for i = 1:numDecompose
    DSall = DSall + DS{i};
end
dsall = find(DSall)';
%	wirebasket
DW = ones(nV, 1);
DW(dsall) = 0;
DW(find(sum(MC(:, dsall), 2))') = 0;
dw = find(DW)'; %	exact wirebasket
%	edges
DE = ones(nV, 1);
DE(dsall) = 0;
DE(dw) = 0;
for i = 1:numDecompose
    DB{i} = double(DB{i} > DW); %	exact boundary list
end
%
DWB = zeros(nV, 1);
DWB(find(sum(MC(:, dw), 2))') = 1;
DWB = double(DWB > DW); %	boundary of wirebasket

%% check the validation of wirebasket
DecomposeUTCheck2(nV, dsall, DE, DW);

%% decomposition information

for i = 1:numDecompose
    numDS(1, i) = sum(DS{i});
end
minSubdomain = min(numDS);
maxSubdomain = max(numDS);
meanSubdomain = floor(sum(numDS) / numDecompose);
numW = sum(DW);
allBoundary = zeros(nV, 1); %	union of boundaries of all decompositions
for i = 1:numDecompose
	allBoundary = allBoundary | DB{i};
end
numallBoundary = sum(allBoundary);
numWB = sum(DWB);

fprintf('decomposition information:\n');
fprintf('number of decomposition: %d \n', numDecompose);
fprintf('minimum size of subdomain in decomposition: %d \n', minSubdomain);
fprintf('maximum size of subdomain in decomposition: %d \n', maxSubdomain);
fprintf('mean size of subdomain in decomposition: %d \n', meanSubdomain);
fprintf('size of all boundaries: %d \n', numallBoundary);
fprintf('size of wirebasket set: %d \n', numW);
fprintf('size of wirebasket set boundary: %d \n\n\n\n', numWB);

%%
end