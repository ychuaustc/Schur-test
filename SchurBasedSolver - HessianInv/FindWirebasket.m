function [DS, DB, DE, DW, DWB] = FindWirebasket(DV, numDecompose, MC, nV, wirebasketType)
%
%   this function construct the wirebasket set for the decomposition in 2 ways:
%   1. take the the edge as the wirebasket set, and the neighbours of the wirebasket set as the the
%   edge
%   2. take all intersection points as the wirebasket set, then add the neighbours of the wirebasket
%   set to the edge 
%
%   INPUT:  DV - decomposition list
%           numDecompose - decomposition number
%           MC - adjacent matrix
%           nV - mesh size
%           wirebasketType - the type of the wirebasket set
%
%   OUTPUT: DS - subdomain list
%           DB - subdomain boundary list
%           DE - edge of all decomposition
%           DW - wirebasket set
%           DWB - wirebasket boundary


fprintf('computing wirebasket set ...\n');
%   start the time recording for the wirebasket computation
tic;


%   pretreatment
for i = 1:numDecompose
    DS{i} = ones(nV, 1);
end
DW = zeros(nV, 1);
DWB = zeros(nV, 1);
DE = zeros(nV, 1);
for i = 1:numDecompose
    for j = 1:numDecompose
        if i ~= j
            dvInd = find(DV{j})';
            DS{i}(dvInd) = 0;
        end
    end
    DB{i} = double(DV{i} > DS{i});
end


%   compute the wirebasket set
switch wirebasketType
    case 1
        for i = 1:numDecompose
            dbInd = find(DB{i})';
            DW(dbInd) = 1;	% get the wirebasket set
        end
        dwInd = find(DW)';
        DWB(find(sum(MC(:, dwInd), 2))') = 1;
        DWB(dwInd) = 0; % get the boundary of the wirebasket set
        DE = DWB;   % get the edge
        DWWB = double(DW | DWB);
        for i = 1:numDecompose
            DS{i} = double(DS{i} > DWWB);   % get the subdomain list
        end
    case 2
        for i = 1:numDecompose
            DE = double(DE | DB{i});
        end
        deInd = find(DE)';
        for i = deInd
            count = 0;
            for j = 1:numDecompose
                if DB{j}(i) ~= 0
                    count = count + 1;
                end
            end
            if count >= 3
                DW(i) = 1;	% get the wirebasket set
            end
        end
        dwInd = find(DW)';
        %
        DWB(find(sum(MC(:, dwInd), 2))') = 1;
        DWB(dwInd) = 0;	% get the boundary of the wirebasket set
        %
        DE = double(DE | DWB);	% get the edge
        DEW = double(DE | DW);
        for i = 1:numDecompose
            DS{i} = double(DS{i} > DEW);    % get the subdomain list
        end
    otherwise
        quit(1)
end


%   finish the time recording for the wirebasket set computation
tw = toc;
fprintf('wirebasket set computed\nrunning time for the computation of the wirebasket: %f s\n\n', tw);


end