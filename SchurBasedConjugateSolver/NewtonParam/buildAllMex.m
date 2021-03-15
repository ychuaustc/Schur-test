mex COMPFLAGS="$COMPFLAGS /openmp /std:c++17" myaccumarray.cpp
mex COMPFLAGS="$COMPFLAGS /openmp" ij2nzIdxs.cpp
mex COMPFLAGS="$COMPFLAGS /MT /std:c++17" -Id:/dev/mkl/include/ -I3rdparty/include/SuiteSparse -I3rdparty/include -L3rdparty/lib/x64 -lSuiteSparse -llinsolvers splsolver_imp.cpp
mex COMPFLAGS="$COMPFLAGS /openmp" -I3rdparty/include projMeshHessians.cpp
mex COMPFLAGS="$COMPFLAGS /openmp" -I3rdparty/include meshIsometricEnergyC.cpp