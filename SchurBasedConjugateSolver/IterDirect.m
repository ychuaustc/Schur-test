function [VertexU, t, iterARAP] = IterDirect(MThetaLVX, MThetaLVY, M, Face, FXV, FYV, FCotTheta, EnULOld, nV, nF, epsArap)
%
%   this function do the iterative step with direct sovler
%
%   INPUT:  MThetaLVX - assign the x-coordinate value of cot(theta) * L(V) to the adjacent matrix
%           MThetaLVY - assign the y-coordinate value of cot(theta) * L(V) to the adjacent matrix
%           M - the coefficient matrix
%           Face - mesh faces
%           FXV - x-coordinates for the face matrix-V
%           FYV - y-coordinates for the face matrix-V
%           FCotTheta - the cotangent value of the face matrix
%           EnULOld - the ARAP energy E(U, L)
%           nV - mesh size
%           nF - number of faces
%           epsArap - convergence error for ARAP method
%
%   OUTPUT: 
%           VertexU - mesh vertices
%           t - running time
%           iterARAP - ARAP method iteration time


iterARAP = 0;   % iteration count

%   iteration
while 1
	iterARAP = iterARAP + 1;
	CX = ArapSystemC(MThetaLVX);
	CY = ArapSystemC(MThetaLVY);
            
    [XU, t_X{iterARAP}] = DirectSolver(M, CX);
    [YU, t_Y{iterARAP}] = DirectSolver(M, CY);
    
	VertexU = [XU, YU];
	[FXU, FYU] = FaceCoord(VertexU, Face);
	[FLVX, FLVY, MThetaLVX, MThetaLVY] ...
    = ArapL(FXV, FYV, FXU, FYU, FCotTheta, Face, nV, nF);	% compute the rotation L for current iteration
 	EnULNew = EnergyUL(FCotTheta, FLVX, FLVY, FXU, FYU);    % compute the ARAP energy for current iteration
    
 	if abs(EnULNew - EnULOld) / nF >= epsArap
        EnULOld = EnULNew;
    else
        break;
	end
end

%   running time
t = 0;
for i = 1:iterARAP
    t = t + t_X{i} + t_Y{i};
end
t = t / (iterARAP * 2);

VertexU = full(VertexU);


end