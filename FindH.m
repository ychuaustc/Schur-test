function [He1, He2, He3, He4] = FindH(x, t)

nv = size(x, 1);
nf = size(t, 1);

faceElens = sqrt( meshFaceEdgeLen2s(x, t) );
faceAngles = meshAnglesFromFaceEdgeLen2(faceElens.^2);
flatXPerFace = [zeros(nf,1) faceElens(:,3) faceElens(:,2).*exp(1i*faceAngles(:,1))];


%% initilization
L = laplacian(x, t, 'uniform');
B = findBoundary(x, t);
I = setdiff(1:nv, B);
Areas = signedAreas(x, t);

isometric_energyies = [ "SymmDirichlet", "ExpSD", "AMIPS", "SARAP", "HOOK", "ARAP", "BARAP", "BCONF"];
hessian_projections = [ "NP", "KP", "FP4", "FP6", "CM" ];

energy_param = 1;
energy_type = isometric_energyies(6);
hession_proj = hessian_projections(2);

findStringC = @(s, names) find(strcmpi(s, names), 1) - 1;
mexoption = struct('energy_type', findStringC(energy_type, isometric_energyies), ...
                   'hessian_projection', findStringC(hession_proj, hessian_projections), ...
                   'energy_param', energy_param, 'verbose', 0);

z = zeros(nv,1);
z(B) = exp(2i*pi*(1:numel(B))'/numel(B));
z(I) = -L(I,I)\(L(I,B)*z(B));

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
[e, g, hs] = meshIsometricEnergyC(fz, gz, D2t, Areas, mexoption);

G = accumarray(g2GIdx(:), g(:));
Hnonzeros = accumarray( nzidx(:), hs(:) ) + Hnonzeros0;
    
%% Newton   
He = replaceNonzeros(H, Hnonzeros);
He1 = He(1:nv, 1:nv);
He2 = He(1:nv, nv + 1:2 * nv);
He3 = He(nv + 1:2 * nv, 1:nv);
He4 = He(nv + 1:2 * nv, nv + 1:2 * nv);
    
end