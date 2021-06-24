function [Hx, Hy, t0, G, dz, x0] = testfunc1(x, xorigin, t, op)

%%
fC2R = @(x) [real(x) imag(x)];
nf = size(t, 1);
nv = size(x, 1);

faceElens = sqrt( meshFaceEdgeLen2s(xorigin, t) );
faceAngles = meshAnglesFromFaceEdgeLen2(faceElens.^2);
flatXPerFace = [zeros(nf,1) faceElens(:,3) faceElens(:,2).*exp(1i*faceAngles(:,1))];


%% initilization
L = laplacian(xorigin, t, 'uniform');
B = findBoundary(xorigin, t);
I = setdiff(1:nv, B);
Areas = signedAreas(xorigin, t);

isometric_energyies = [ "SymmDirichlet", "ExpSD", "AMIPS", "SARAP", "HOOK", "ARAP", "BARAP", "BCONF"];
hessian_projections = [ "NP", "KP", "FP4", "FP6", "CM" ];

energy_param = 1;
energy_type = isometric_energyies(1);
hession_proj = hessian_projections(2);

findStringC = @(s, names) find(strcmpi(s, names), 1) - 1;
mexoption = struct('energy_type', findStringC(energy_type, isometric_energyies), ...
                   'hessian_projection', findStringC(hession_proj, hessian_projections), ...
                   'energy_param', energy_param, 'verbose', 0);

z = zeros(nv,1);
if op == 1
    z(B) = exp(2i*pi*(1:numel(B))'/numel(B));
    z(I) = -L(I,I)\(L(I,B)*z(B));
else
    z = x(:, 1) + 1i * x(:, 2);
end
x0(:, 1:2) = fC2R(z);
x0(:, 3) = 0;

figure; h = trimesh(t, real(z), imag(z), z*0); hold on; view(2); axis equal; axis off;
set(h, 'FaceColor', 'w', 'edgealpha', 0.1, 'edgecolor', 'k');

D2 = -1i/4*(flatXPerFace(:,[2 3 1])-flatXPerFace(:,[3 1 2]))./Areas;
D = sparse(repmat(1:nf,3,1)', t, D2);
D2t = D2.';

fDeformEnergy = @(z) meshIsometricEnergyC(conj(D*conj(z)), D*z, D2t, Areas, mexoption); 

en = fDeformEnergy(z);

lambda = 1e-8;
%% initialization, get sparse matrix pattern
[xmesh, ymesh] = meshgrid(1:6, 1:6);
t2 = [t t+nv]';
Mi = t2(xmesh(:), :);
Mj = t2(ymesh(:), :);

% H = sparse(Mi, Mj, 1, nv*2, nv*2);  % only pattern is needed
H = [L L; L L];


% nonzero indices of the matrix
Hnonzeros0 = zeros(nnz(H),1);
idxDiagH = ij2nzIdxs(H, uint64(1:nv*2), uint64(1:nv*2));
Hnonzeros0(idxDiagH) = lambda*2;
nzidx = ij2nzIdxs(H, uint64(Mi), uint64(Mj));


%% main loop
g2GIdx = uint64(t2);

fz = conj(D*conj(z)); % equivalent but faster than conj(D)*z;
gz = D*z;
[~, g, hs] = meshIsometricEnergyC(fz, gz, D2t, Areas, mexoption);

G = accumarray(g2GIdx(:), g(:));
Hnonzeros = accumarray( nzidx(:), hs(:) ) + Hnonzeros0;
    
%% Newton
% H = sparse(Mi, Mj, hs, nv*2, nv*2) + 2*lambda*sparse(1:nv*2, 1:nv*2, 1, nv*2, nv*2);    
H = replaceNonzeros(H, Hnonzeros);
Hx = H(1:nv, 1:nv);
Hy = H(nv + 1:2 * nv, nv + 1:2 * nv);
tic;
dz0 = H \ -G;
t0 = toc;
dz(1:nv) = Hx \ -G(1:nv);
dz(nv + 1:2 * nv) = Hy \ -G(nv + 1:2 * nv);

end