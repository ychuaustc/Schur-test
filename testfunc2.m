function [iter, t, dz] = testfunc2(MSSx, MSEx, MSEallx, MWWx, MWEx, MEEx, MSSy, MSEy, MSEally, MWWy, MWEy, MEEy, ...
                                   dsInd, dwInd, deInd, dssize, vwe, vweorigin, fwe, ide, G, ...
                                   nV, numDecompose, epsSchur)

LW = 1;
LE = 1;

[Hxd, Hyd, ~, ~, ~] = testfunc1(vwe, vweorigin, fwe, 2);
MEE_Wx = Hxd(ide, ide);
MEE_Wy = Hyd(ide, ide);

MSSsolverx = findsolverSS(MSSx, numDecompose);
MSSsolvery = findsolverSS(MSSy, numDecompose);
presolverx = findsolverPre(MWWx, MWEx, MEE_Wx);
presolvery = findsolverPre(MWWy, MWEy, MEE_Wy);

CX = -G(1:nV);
CY = -G(nV + 1:2 * nV);
[CSX, CWX, bX] = SchurSystemC(MSEallx, MWWx, MWEx, CX, dsInd, dwInd, deInd, MSSsolverx, numDecompose);
[CSY, CWY, bY] = SchurSystemC(MSEally, MWWy, MWEy, CY, dsInd, dwInd, deInd, MSSsolvery, numDecompose);

% [XU, iter_X, t_X] ...
% = SchurConjSolver(MSSx, MSEx, MWWx, MWEx, MEEx, CSX, CWX, bX, dsInd, dwInd, deInd, nV, numDecompose, epsSchur);
% [YU, iter_Y, t_Y] ...
% = SchurConjSolver(MSSy, MSEy, MWWy, MWEy, MEEy, CSY, CWY, bY, dsInd, dwInd, deInd, nV, numDecompose, epsSchur);
[XU, iter_X, t_X] ...
= SchurConjPreSolver(MSEallx, MWWx, MWEx, MEEx, MEE_Wx, LW, LE, CSX, CWX, bX, ...
                     dsInd, dwInd, deInd, MSSsolverx, presolverx, dssize, nV, numDecompose, epsSchur);
[YU, iter_Y, t_Y] ...
= SchurConjPreSolver(MSEally, MWWy, MWEy, MEEy, MEE_Wy, LW, LE, CSY, CWY, bY, ...
                     dsInd, dwInd, deInd, MSSsolvery, presolvery, dssize, nV, numDecompose, epsSchur);
dz(1:nV) = XU;
dz(nV + 1:2 * nV) = YU;

t = t_X + t_Y;
iter = iter_X + iter_Y;