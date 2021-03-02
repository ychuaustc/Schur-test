function [DS, DE, DW] = Decomp(MC, Vertex, nV, numDecompose)
%
%   this function decompose the mesh into #numDecompose pieces, each piece is obtained in
%   an iterative process
%
%   INPUT:  MC - adjacent matrix
%           Vertex - mesh vertices
%           nV - mesh size
%           numDecompose - decomposition number
%
%   OUTPUT: DS - subdomain list
%           DE - edge of all decomposition
%           DW - wirebasket set


%   decompose the mesh using dice
Vertex = Vertex(:, 1:2);
level = log2(numDecompose);
[Map, sepEdge] = dice('geopart', level, MC, Vertex, 30);

%   plot the decomposed mesh
gplotmap_copy(MC, Vertex, Map);

%   find the subdomains and edge / wirebasket
%   subdomains
parts = unique(Map);
edge = unique(sepEdge);
for itp = parts
    part = find(Map == itp)';
    dv = zeros(nV, 1);
    dv(part) = 1;
    dv(edge) = 0;
    DS{itp + 1} = dv;
end
%   edge / wirebasket
subdomain = setdiff(1:nV, edge)';
w = setdiff(1:nV, find(sum(MC(:, subdomain), 2)));
w = setdiff(w, subdomain);  % or w = intersect(w, edge);
DW = zeros(nV, 1);
DW(w) = 1;
DE = zeros(nV, 1);
DE(edge) = 1;
DE(w) = 0;

%   display the details of the decomposition
for i = 1:numDecompose
    numDS(1, i) = sum(DS{i});
end
minSubdomain = min(numDS);
maxSubdomain = max(numDS);
meanSubdomain = floor(sum(numDS) / numDecompose);
numW = sum(DW);
numE = sum(DE);

fprintf('decomposition informations:\n');
fprintf('the number of decompositions: %d \n', numDecompose);
fprintf('the minimum size of all subdomains: %d \n', minSubdomain);
fprintf('the maximum size of all subdomains: %d \n', maxSubdomain);
fprintf('the mean size of all subdomains: %d \n', meanSubdomain);
fprintf('the size of the edge and wirebasket: %d \n', numE + numW);
fprintf('mesh size / edge-wirebasket size: %f \n\n\n', nV / (numE + numW));


end