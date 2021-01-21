function [U11, U12, U21, U22, V11, V12, V21, V22] = SVD2D(A, B, C, D)
%
%   this function computes the SVD of (A B;C D)
%
%   INPUT:  
%
%   OUTPUT: 


    E = (A + D) / 2;
    F = (A - D) / 2;
    G = (C + B) / 2;
    H = (C - B) / 2;
    Q = sqrt(E.^2 + H.^2);
    R = sqrt(F.^2 + G.^2);
    S(:, 1) = Q + R;
    S(:, 2) = Q - R;
    a1 = atan2(G, F);
    a2 = atan2(H, E);
    theta = (a2 - a1) / 2;
    phi = (a2 + a1) / 2;
    U11 = cos(phi);
    U21 = sin(phi);
    U12 = -sin(phi);
    U22 = cos(phi);
    V11 = cos(theta);
    V21 = -sin(theta);
    V12 = sin(theta);
    V22 = cos(theta);
    
    
end