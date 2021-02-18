function [DS, DE, DW, DWB] = Decomp(MC, Vertex, nV, numDecompose)
%
%   this function decompose the mesh into #numDecompose pieces, each piece is obtained in
%   an iterative process
%
%   INPUT:  MC - adjacent matrix
%           Vertex - mesh vertices
%           numDecompose - decomposition number
%
%   OUTPUT: DV - decomposition list
%           DS - subdomain list
%           DB - subdomain boundary list
%           DE - edge of all decomposition
%           DW - wirebasket set
%           DWB - wirebasket boundary


%   decompose the mesh using dice
Vertex = Vertex(:, 1:2);
level = log2(numDecompose);
[Map, sepEdge] = dice('geopart', level, MC, Vertex, 30);

%   plot the decomposed mesh
gplotmap(MC, Vertex, Map);

%   find the subdomains and edge / wirebasket
%   subdomains
parts = unique(Map);
edge = unique(sepEdge);
for itp = parts
    part = find(Map == itp)';
    dv = zeros(nV, 1);
    dv(part) = 1;
    DS{itp + 1} = dv;
    DS{itp + 1}(edge) = 0;
end
%   edge / wirebasket
subdomain = setdiff(1:nV, edge)';
w = setdiff(1:nV, find(sum(MC(:, subdomain), 2)));
DW = zeros(nV, 1);
DW(w) = 1;
DE = zeros(nV, 1);
DE(edge) = 1;
DE(w) = 0;
DWB = DE;
for itp = parts
    DS{itp + 1}(w) = 0;
end


end