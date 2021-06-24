function [DS, DE, DV, numDS, numDecompose] = Decompose(nV, nv, MC)
%
%   this function decompose the mesh into #numDecompose pieces, each piece is obtained in
%   an iterative process
%
%   INPUT:  nV - mesh size
%           nv - expected size of each decomposition
%           B - boundary
%           MC - adjacent matrix
%           Vertex - mesh vertices
%           Face - mesh faces
%           nF - number of faces
%
%   OUTPUT: DV - decomposition list
%           DS - subdomain list
%           DB - subdomain boundary list
%           DE - edge of all decomposition
%           DW - wirebasket set
%           DWB - wirebasket boundary
%           numDS - subdomain size list
%           numDecompose - decomposition number


%%	initialization
ML = MC;	% adjacent matrix for the left part after one decomposition step
verticeL = ones(nV, 1);	% vertice in the left part after one decomposition step
column = 1;	% step count


%%	decomposition loops
while sum(verticeL) > 0
    %	pick up a random vertex on B, add it to vertice
    boundTemp = find(boundL);
    verticeAdd = zeros(nV, 1);
    vertice = zeros(nV, 1);
    verticeAdd(boundTemp(1)) = 1;
    vertice(boundTemp(1)) = 1;
    num = 1;
    
    %   find one decomposition iteratively, and save the result to vertice
    [~, vertice, ~, ~] = DecomposeIter(verticeAdd, vertice, ML, num, nV, nv);
    
    %	compute the adjacent matrix for this decomposition
    verticeLWOI = verticeL - vertice;   % WOI means without interface
    verticeLWOITemp = find(verticeLWOI)';
    MD = ML;	% adjacent matrix for the decomposition
    MD(verticeLWOITemp, :) = 0;
    MD(:, verticeLWOITemp) = 0;
    
    %   compute the interface of the decomposition and the left part
    MLWOI = ML - MD;
    verticeTemp = find(vertice)';
    InterfaceDTemp = verticeTemp(find(sum(MLWOI(:, verticeTemp), 1))');
    
    %	modification of the decomposition, the adjacent matrix and the interface: if a vertex is in the left part, and 
    %	is only connected to the decomposition, we add it to the decomposition
    interfaceDNeighborTemp = find(sum(ML(:, InterfaceDTemp), 2))';
    interfaceDNeighborComp = ones(nV, 1);
    interfaceDNeighborComp(InterfaceDTemp) = 0;
    interfaceDNeighborCompTemp = find(interfaceDNeighborComp)';
    NeighborToAddTemp = interfaceDNeighborTemp(find(sum(ML(interfaceDNeighborCompTemp, interfaceDNeighborTemp), 1) == 0)');
    NeighborToAdd = zeros(nV, 1);
    NeighborToAdd(NeighborToAddTemp) = 1;
    NeighborToAddTemp = find(double(NeighborToAdd > vertice))';
    vertice(NeighborToAddTemp) = 1;	% here we get the "FINAL" decomposition set IN THE ITERATIVE PROCESS
    
    %   compute the adjacenet matrix
    verticeLWOI = verticeL - vertice;
    verticeLWOITemp = find(verticeLWOI)';
    MD = ML;
    MD(verticeLWOITemp, :) = 0;
    MD(:, verticeLWOITemp) = 0;	% here we get the "FINAL" adgjacent matrix IN THE ITERATIVE PROCESS
    
    %   compute the interface of the decomposition
    MLWOI = ML - MD;
    verticeTemp = find(vertice)';
    InterfaceDTemp = verticeTemp(find(sum(MLWOI(:, verticeTemp), 1))');
    InterfaceD = zeros(nV, 1);
    InterfaceD(InterfaceDTemp) = 1; % here we get the "FINAL" interface (of the decomposition and the left part) 
                                    % IN THE ITERATIVE PROCESS
    boundD = InterfaceD | (boundL & vertice); %	here we get the "FINAL" boundary IN THE ITERATIVE PROCESS
    
    %	construct the decomposition list
    DV{column} = vertice; %	the decomposition
    DB{column} = boundD;  %	the boundary of the decomposition
    DM{column} = MD;      %	the adjacent matrix of the decomposition
    
    %	compute the boundary, adjacent matrix and vertice for the left part
    boundL = (boundL > vertice) | InterfaceD;	% the boundary of the left part
    MI = ML;
    verticeWOI = verticeL - InterfaceD;
    verticeWOITemp = find(verticeWOI)';
    MI(verticeWOITemp, :) = 0;
    MI(:, verticeWOITemp) = 0;
    ML = MLWOI + MI;	% the adjacent matrix of the left part
    verticeLWOI = verticeL - vertice;
    verticeL = verticeLWOI | InterfaceD;	% the vertice of the left part
    
    %	update the count of the loop
    column = column + 1;
end


%%	merge small decomposition into large ones
numDecompose = size(DV, 2);
nMerge = ceil(nv / 5);	% if the size of one decomposition does not reach this limit, we merge it into another one
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
numDecompose = size(DV, 2);


%%  finish time recording for decomposition
td = toc;


%%	check the validation of decomposition
DecomposeCheck(nV, B, DV, DB, numDecompose, td);


%%  construct the wirebasket set
[DS, DB, DE, DW, DWB] = FindWirebasket(DV, numDecompose, MC, nV);


%% display the details of the decomposition
for i = 1:numDecompose
    numDS(1, i) = sum(DS{i});
end
minSubdomain = min(numDS);
maxSubdomain = max(numDS);
meanSubdomain = floor(sum(numDS) / numDecompose);
numW = sum(DW);
numE = sum(DE);
numWB = sum(DWB);

fprintf('decomposition informations:\n');
fprintf('the number of decompositions: %d \n', numDecompose);
fprintf('the minimum size of all subdomains: %d \n', minSubdomain);
fprintf('the maximum size of all subdomains: %d \n', maxSubdomain);
fprintf('the mean size of all subdomains: %d \n', meanSubdomain);
fprintf('the size of the edge (the union of the boundaries of all subdomains): %d \n', numE);
fprintf('the size of the wirebasket set: %d \n', numW);
fprintf('the size of the boundary of the wirebasket: %d \n\n\n', numWB);


%%  draw the decomposed mesh
% DrawDecomposedMesh(Vertex, Face, DV, nV, nF, numDecompose);


end