function solverMSS = SolverCompute(MSS, numDecompose)

for i = 1:numDecompose
    solverMSS{i} = splsolver(MSS{i}, 'ldlt');
end