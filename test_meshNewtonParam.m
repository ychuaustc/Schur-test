clear;
clc;
profile on;
%
fC2R = @(x) [real(x) imag(x)];
fR2C = @(x) complex(x(:,1), x(:,2));
%
PathSet();
%
[meshType, nV, nv, numDecompose, fileName, epsArap, epsSchur, NewtonOp] = SetParameter();
[Vertex, Face] = MeshGeneration(meshType, nV, fileName);
[nF, I, nI, B, nB, MC] = MeshInfo(Vertex, Face, nV);
%
[DS, DE, DW, Map, sepEdge] = Decomp(MC, Vertex, nV, numDecompose);
[dsInd, dwInd, deInd, dsIndall, dweInd, ide, dssize] = SchurSystemIndex(DS, DE, DW, numDecompose);
%
vweorigin = Vertex(dweInd, :);
[FinWE, finweInd] = FindFinWE(Face, dweInd, nV);
fwe = FindFforWE(FinWE, dweInd);


%%
faceElens = sqrt( meshFaceEdgeLen2s(Vertex, Face) );
faceAngles = meshAnglesFromFaceEdgeLen2(faceElens.^2);
flatXPerFace = [zeros(nF,1) faceElens(:,3) faceElens(:,2).*exp(1i*faceAngles(:,1))];


%% initilization
L = laplacian(Vertex, Face, 'uniform');
Areas = signedAreas(Vertex, Face);

isometric_energyies = [ "SymmDirichlet", "ExpSD", "AMIPS", "SARAP", "HOOK", "ARAP", "BARAP", "BCONF"];
hessian_projections = [ "NP", "KP", "FP4", "FP6", "CM" ];

energy_param = 1;
energy_type = isometric_energyies(1);
hession_proj = hessian_projections(2);

findStringC = @(s, names) find(strcmpi(s, names), 1) - 1;
mexoption = struct('energy_type', findStringC(energy_type, isometric_energyies), ...
                   'hessian_projection', findStringC(hession_proj, hessian_projections), ...
                   'energy_param', energy_param, 'verbose', 0);

z = zeros(nV,1);
z(B) = exp(2i*pi*(1:numel(B))'/numel(B));
z(I) = -L(I,I)\(L(I,B)*z(B));

figure; h = trimesh(Face, real(z), imag(z), z*0); hold on; view(2); axis equal; axis off;
set(h, 'FaceColor', 'w', 'edgealpha', 0.1, 'edgecolor', 'k');

D2 = -1i/4*(flatXPerFace(:,[2 3 1])-flatXPerFace(:,[3 1 2]))./Areas;
D = sparse(repmat(1:nF,3,1)', Face, D2);
D2t = D2.';

fDeformEnergy = @(z) meshIsometricEnergyC(conj(D*conj(z)), D*z, D2t, Areas, mexoption); 

en = fDeformEnergy(z);

lambda = 1e-8;
%% initialization, get sparse matrix pattern
[xmesh, ymesh] = meshgrid(1:6, 1:6);
t2 = [Face Face+nV]';
Mi = t2(xmesh(:), :);
Mj = t2(ymesh(:), :);

H = [L L; L L];

% nonzero indices of the matrix
Hnonzeros0 = zeros(nnz(H),1);
idxDiagH = ij2nzIdxs(H, uint64(1:nV*2), uint64(1:nV*2));
Hnonzeros0(idxDiagH) = lambda*2;
nzidx = ij2nzIdxs(H, uint64(Mi), uint64(Mj));


%% main loop
g2GIdx = uint64(t2);
for it=1:1000
    tt = tic;

    fz = conj(D*conj(z)); % equivalent but faster than conj(D)*z;
    gz = D*z;
    [e, g, hs] = meshIsometricEnergyC(fz, gz, D2t, Areas, mexoption);

    G = accumarray(g2GIdx(:), g(:));
    Hnonzeros = accumarray( nzidx(:), hs(:) ) + Hnonzeros0;
    
    %% Newton
    He = replaceNonzeros(H, Hnonzeros);
    
    switch NewtonOp
        case 1
            dz = He \ (-G);
        case 2
            Hx = He(1:nV, 1:nV);
            Hy = He(nV + 1:2 * nV, nV + 1:2 * nV);
            dz(1:nV) = Hx \ -G(1:nV);
            dz(nV + 1:2 * nV) = Hy \ -G(nV + 1:2 * nV);
        case 3
            LW = 1;
            LE = 1;
            Hx = He(1:nV, 1:nV);
            Hy = He(nV + 1:2 * nV, nV + 1:2 * nV);
            [MSSx, MSEx, MSEallx, MWWx, MWEx, MEEx] = SchurSystemM(Hx, dsInd, dwInd, deInd, dsIndall, numDecompose);
            [MSSy, MSEy, MSEally, MWWy, MWEy, MEEy] = SchurSystemM(Hy, dsInd, dwInd, deInd, dsIndall, numDecompose);
            %
            Vertex0(:, 1:2) = fC2R(z);
            Vertex0(:, 3) = 0;
            vwe = Vertex0(dweInd, :);
            %
            [Hxd, Hyd] = FindHessian(vwe, vweorigin, fwe);
            MEE_Wx = Hxd(ide, ide);
            MEE_Wy = Hyd(ide, ide);
            MSSsolverx = findsolverSS(MSSx, numDecompose);
            MSSsolvery = findsolverSS(MSSy, numDecompose);
%             MWWsolverx = findsolverW(MWWx);
%             MWWsolvery = findsolverW(MWWy);
            presolverx = findsolverPre(MWWx, MWEx, MEE_Wx);
            presolvery = findsolverPre(MWWy, MWEy, MEE_Wy);
            CX = -G(1:nV);
            CY = -G(nV + 1:2 * nV);
            [CSX, CWX, bX] = SchurSystemC(MSEallx, MWWx, MWEx, CX, dsInd, dwInd, deInd, MSSsolverx, numDecompose);
            [CSY, CWY, bY] = SchurSystemC(MSEally, MWWy, MWEy, CY, dsInd, dwInd, deInd, MSSsolvery, numDecompose);
            [XU, iter_X{it}, t_X{it}] ...
            = SchurConjPreSolver(MSEallx, MWWx, MWEx, MEEx, MEE_Wx, LW, LE, CSX, CWX, bX, ...
                                 dsInd, dwInd, deInd, MSSsolverx, presolverx, dssize, nV, numDecompose, epsSchur);
            [YU, iter_Y{it}, t_Y{it}] ...
            = SchurConjPreSolver(MSEally, MWWy, MWEy, MEEy, MEE_Wy, LW, LE, CSY, CWY, bY, ...
                                 dsInd, dwInd, deInd, MSSsolvery, presolvery, dssize, nV, numDecompose, epsSchur);
            dz(1:nV) = XU;
            dz(nV + 1:2 * nV) = YU;
        otherwise
            quit(1);
    end

    dz = fR2C( reshape(dz, [], 2) );
    
    %% orientation preservation
    ls_t = min( maxtForPositiveArea( fz, gz, conj(D*conj(dz)), D*dz )*0.9, 1 );

    %% line search energy decreasing
    fMyFun = @(t) fDeformEnergy( dz*t + z );
    normdz = norm(dz);

    dgdotfz = dot( G, [real(dz); imag(dz)] );
    
    ls_alpha = 0.2; ls_beta = 0.5;
    fQPEstim = @(t) en+ls_alpha*t*dgdotfz;

    e_new = fMyFun(ls_t);
    while ls_t*normdz>1e-6 && e_new > fQPEstim(ls_t)
        ls_t = ls_t*ls_beta;
        e_new = fMyFun(ls_t);
    end
    if ls_t*normdz<=1e-1
        break;
    end
    en = e_new;
    
    tt0{it} = toc(tt);
    fprintf('it: %3d, en: %.3e, runtime: %fs, ls: %.2e\n', it, en, tt0{it}, ls_t);
    
    %% update
    z = dz*ls_t + z;

    %%
    title( sprintf('iter %d', it) );
    set(h, 'Vertices', fC2R(z));
    drawnow;
    pause(0.001);
end

%%
profile viewer;