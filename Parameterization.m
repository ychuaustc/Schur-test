%%  return parameterized mesh
%%
%%  Input:  //
%%  Output: //

function [vertex] = Parameterization(M, Cx, Cy, Vertex, I, B, nB)

%% do the parameterization by solving the linear system
fprintf('starting parameterization ...\n');
tic;
X = (M \ Cx);
Y = (M \ Cy);
fprintf('parameterization done\nparameterization running time: %f s \n\n\n\n', toc / 2);

%% set the parameterized meshes
vertex = Vertex;
vertex(:, 3) = 0;
vertex(B, 1:2) = SquareBoundary(nB);
vertex(I, 1) = X;
vertex(I, 2) = Y;

end