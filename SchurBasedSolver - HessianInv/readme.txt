This is a test that checks whether the Schur complement method works for the second order approximation method. 
For the Hessian matrix H of the energy of the deformation, we simply apply the Schur method to compute H^-1 V 
for a given V in R(H), and see how much can we reduce the iterations in the Schur complement method with a 
preconditioner.

to use this code,

-run HessianInv.m
-select the parameters
-see the results


files included:

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


This is just a prototype for test and more modification is needed to see how this method works, such as test on more 
type of energies, parallelization, and so on. For time being, only ARAP, symmetric Dirichlet and SARAP are tested 
for some meshes. For the tests we've done, the iteration time can be 1 / 6 compared to the un-preconditioned case 
for most cases.