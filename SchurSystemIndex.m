function [dsInd, dwInd, deInd, dsIndall, dweInd, ide, dssize] = SchurSystemIndex(DS, DE, DW, numDecompose)
%
%   this function computes the index and matrix blocks for the Schur system upon a given
%   decomposition
%
%   INPUT:  DS - subdomain list
%           DE - edge of all decompositions
%           DW - wirebasket set
%           numDecompose - decomposition number
%
%   OUTPUT: dsInd - the index for the subdomains
%           dwInd - the index for the wirebasket
%           deInd - the index for the edge


%   the index in the decomposition
for i = 1:numDecompose
    dsInd{i} = find(DS{i});	% the index for the subdomains
end
dwInd = find(DW);	% the index for the wirebasket
deInd = find(DE);	% the index for the edge
dweInd = [dwInd;deInd];

dsIndall = [];
k = 1;
for i = 1:numDecompose
    dsIndall = [dsIndall; dsInd{i}];
    dssize{i} = [k: k + size(dsInd{i}, 1) - 1]';
    k = k + size(dsInd{i});
end

nw = size(dwInd, 1);
ne = size(deInd, 1);
ide = [nw + 1:nw + ne]';

end