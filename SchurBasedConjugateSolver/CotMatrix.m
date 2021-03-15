function [FCotTheta, MCotTheta] = CotMatrix(Face, FX, FY, nV)
%
%   this function computes the cotangent matrix for the mesh, which will be used to compute the ARAP
%   system
%
%   INPUT:  Face - mesh faces
%           FX - x-coordinates for the face matrix
%           FY - y-coordinates for the face matrix
%           nV - mesh size
%
%   OUTPUT: FCotTheta - the cotangent value of the face matrix
%           MCotTheta - the cotangent value of the adjacent matrix


F1X = FX(:, 1);
F2X = FX(:, 2);
F3X = FX(:, 3);
F1Y = FY(:, 1);
F2Y = FY(:, 2);
F3Y = FY(:, 3);

E1X = F3X - F2X;
E1Y = F3Y - F2Y;
E2X = F1X - F3X;
E2Y = F1Y - F3Y;
E3X = F2X - F1X;
E3Y = F2Y - F1Y;
lE1 = sqrt(E1X .^ 2 + E1Y .^ 2);	% length of the edge toward vertex 1 in triangle
lE2 = sqrt(E2X .^ 2 + E2Y .^ 2);    % length of the edge toward vertex 2 in triangle
lE3 = sqrt(E3X .^ 2 + E3Y .^ 2);    % length of the edge toward vertex 3 in triangle

sinTheta1 = abs(E2X .* E3Y - E3X .* E2Y) ./ (lE3 .* lE2);
sinTheta2 = abs(E3X .* E1Y - E1X .* E3Y) ./ (lE1 .* lE3);
sinTheta3 = abs(E1X .* E2Y - E2X .* E1Y) ./ (lE2 .* lE1);
cosTheta1 = -(E3X .* E2X + E3Y .* E2Y) ./ (lE3 .* lE2);
cosTheta2 = -(E1X .* E3X + E1Y .* E3Y) ./ (lE1 .* lE3);
cosTheta3 = -(E2X .* E1X + E2Y .* E1Y) ./ (lE2 .* lE1);
cotTheta1 = cosTheta1 ./ sinTheta1;
cotTheta2 = cosTheta2 ./ sinTheta2;
cotTheta3 = cosTheta3 ./ sinTheta3;
FCotTheta = [cotTheta1 cotTheta2 cotTheta3];

MCotTheta = sparse(Face, Face(:, [2 3 1]), [cotTheta3 cotTheta1 cotTheta2], nV, nV);


end