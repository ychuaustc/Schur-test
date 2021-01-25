This is a test that check whether the Schur complement method works for the Hessian matrix, it simply computes the Hessian of the ARAP energy, say H, then apply the Schur method to compute H^-1 V for any V in R(H), and see whether the preconditioned can help reduce the iteration time. As shown in the result, the reduction can be 1 / 6 for unstructured mesh. 

run:

-run Parameterization.m
-select the parameters
-see the results


files:

-HessianInv.m
-SetParameter.m
-MeshGeneration.m
      -GenerateST.m
      -GenerateUT.m
      -GenerateFromObj.m
-MeshInfo.m
-Decompose.m
      -DecomposeIter.m
      -FindWirebasket.m
      -DecomposeCheck.m
-getM.m
-SchurSystemM.m
-getH.m
-DirectSolver.m
-SchurSystemC.m
-SchurConjSolver.m
-SchurConjPreSolver.m
      -SchurMultiply.m
      -Preconditioner.m


This is just a prototype for test and more modification is needed to see how this method works, such as test on other energies, parallelization, and so on.