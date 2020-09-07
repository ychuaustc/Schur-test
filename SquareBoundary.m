%%  return consequent coordinates on unit square for given boundary
%%
%%  Input:  nB: boundary size (constant)
%%  Output: coord: 2d coordinates consequantly located on unit square, start from (0, 0) (nB * 2 matrix)

function coord = SquareBoundary(nB)

nPerSide = ceil(nB / 4);
nLastSide = nB - nPerSide * 3;
    
side = linspace(0, 1, nPerSide + 1);
lastside = linspace(0, 1, nLastSide + 1);
    
side1 = zeros(nPerSide, 2);
side2 = zeros(nPerSide, 2);
side3 = zeros(nPerSide, 2);
side4 = zeros(nLastSide, 2);
    
side1(:, 1) = side(1:nPerSide);
side1(:, 2) = 0;
    
side2(:, 1) = 1;
side2(:, 2) = side(1:nPerSide);
    
side3(:, 1) = side(nPerSide + 1:-1:2);
side3(:, 2) = 1;
    
side4(:, 1) = 0;
side4(:, 2) = lastside(nLastSide + 1:-1:2);
    
coord = [side1;side2;side3;side4];

end