run:

-run Parameterization.m
-select the parameters
-see the results


files:

-Parameterization.m
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
      -DrawDecomposedMesh.m
-FaceCoord.m
-CotMatrix.m
-Initialization.m
-ArapSystemM.m
-SchurSystemM.m
-LagMultip.m
-MWirebasket.m
-ArapL.m
      -SVD2D.m
-EnergyUL.m
-ArapSystemC.m
-DirectSolver.m
-SchurSystemC.m
-SchurConjSolver.m
-SchurConjPreSolver.m
      -SchurMultiply.m
      -Preconditioner.m
-DrawMesh.m


note

-


to be added / modified / tested ...

-the decomposition algorithm has to be updated
-take wirebaskt set as the neighbor of all subdomain intersections
-more examples is awaiting to be tested (ball, face, bunny, etc)
-find a good initialization
-optimization is needed
-parallelization is needed
-takes too long to draw the decomposed mesh