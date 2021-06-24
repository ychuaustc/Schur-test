#define _CRT_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS

#include <mex.h>
#include <string>
#include <vector>

#include "splsolver.h"

const mxClassID mxPOINTER_CLASS = mxINDEX_CLASS;
using pointer_t = size_t;


//////////////////////////////////////////////////////////////////////////
/* memory management - redirect to MATLAB memory manager */
void* operator new(size_t size)
{
    void* ptr = mxMalloc(size);
    mexMakeMemoryPersistent(ptr);
    return ptr;
}
void* operator new[](size_t size)
{
    void* ptr = mxMalloc(size);
    mexMakeMemoryPersistent(ptr);
    return ptr;
}

void operator delete(void* ptr) { mxFree(ptr); }
void operator delete[](void* ptr) { mxFree(ptr); }


//////////////////////////////////////////////////////////////////////////
// replace mxGetNzmax with nnz, as mxGetNzmax has over-allocated elements (zeros)
size_t nnz(const mxArray* m) { return *(mxGetJc(m) + mxGetN(m)); }

void verifyMatIsSquareSparse(const mxArray* m, bool checkDoublePrecision = true)
{
    if (!mxIsSparse(m)) mexErrMsgTxt("Invalid input: matrix should be sparse.");
    if (mxGetM(m) != mxGetN(m)) mexErrMsgTxt("Invalid input: matrix should be square matrix (nxn).");

    if (checkDoublePrecision && !mxIsDouble(m)) mexErrMsgTxt("Invalid input: matrix should be in double precision.");
}

Eigen::Map<const SpMat> extractSpMat(const mxArray* m, double* replaceNonzeros = nullptr)
{
    static_assert(sizeof(index_t) == 8, "only for 64bit");
    return Eigen::Map<const SpMat>(mxGetM(m), mxGetN(m), nnz(m), (const index_t*)mxGetJc(m), (const index_t*)mxGetIr(m), replaceNonzeros ? replaceNonzeros : mxGetPr(m));
}

enum spl_solve_stage {
    LST_SYMFACT = 1,
    LST_NUMFACT = 2,
    LST_SOLVE = 3,
    LST_NUMFACT_AND_SOLVE = 4,
    LST_ALL = 5,
    LST_CLOSE_HANDLE = 6
};

//////////////////////////////////////////////////////////////////////////
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    if (nrhs < 2)
        mexErrMsgTxt("Invalid input: not enough input, x = SPLsolver(mode, handle, A, Av, b);");

    using MapMatC = Eigen::Map<const Eigen::MatrixXd>;
    using MapMat = Eigen::Map<Eigen::MatrixXd>;

    const int mode = int(mxGetScalar(prhs[0]));

    if (mode == LST_SYMFACT) {
        // SPLsolver(1, A, solvertype)
        //check input var
        const mxArray* cell_A = prhs[1];
        if (!mxIsCell(cell_A)) mexErrMsgTxt("Invalid input: A must be a cell!");
        for (int i = 0; i < mxGetNumberOfElements(cell_A); i++)
        {
            const mxArray* mat = mxGetCell(cell_A, i);
            verifyMatIsSquareSparse(mat); //check precision = double
        }
        //create & return handle
        std::vector<SparseLinearSolver*>* solver = new std::vector<SparseLinearSolver*>;
        plhs[0] = mxCreateNumericMatrix(1, 1, mxPOINTER_CLASS, mxREAL);
        pointer_t* ptr_handle = (pointer_t*)mxGetData(plhs[0]);
        *ptr_handle = (pointer_t)solver;
 
        //create solver & analyze pattern
        int solvertype = nrhs > 2 ? int(mxGetScalar(prhs[2])) : SLS_LU;  //default mode = LU    
        solver->resize(mxGetNumberOfElements(cell_A));
#pragma omp parallel for
        for (int i = 0; i < mxGetNumberOfElements(cell_A); i++)
        {
            solver->at(i) = new SparseLinearSolver(solvertype);
            const mxArray* mat = mxGetCell(cell_A, i);
            solver->at(i)->analyzePattern(extractSpMat(mat));
        }
        return;
    }

    if (mode == LST_ALL) {
        // SPLsolver(5, A, b, solvername)
        //check input var
        if (nrhs < 3)
            mexErrMsgTxt("Invalid input: not enough input, SPLsolver(mode, A, b, solvername)");

        const mxArray* cell_A = prhs[1];
        const mxArray* cell_b = prhs[2];
        if (!mxIsCell(cell_A) || !mxIsCell(cell_b))
            mexErrMsgTxt("Invalid input: A and b must be cells!");
        if (mxGetNumberOfElements(cell_A) != mxGetNumberOfElements(cell_b))
            mexErrMsgTxt("Invalid input: size(A) must equal to size(b)!");
        for (int i = 0; i < mxGetNumberOfElements(cell_A); i++)
        {
            const mxArray* mat = mxGetCell(cell_A, i);
            verifyMatIsSquareSparse(mat); //check precision = double
            const mxArray* vec = mxGetCell(cell_b, i);
            if (!mxIsDouble(vec)) mexErrMsgTxt("Invalid input: vector b should be in double precision.");
            if (mxGetM(vec) != mxGetM(mat)) mexErrMsgTxt("Invalid input: the size of vector b does not match the coefficient matrix.");
        }

        //solve
        plhs[0] = mxDuplicateArray(cell_b);
        if (nlhs > 1) plhs[1] = mxCreateLogicalMatrix(mxGetNumberOfElements(cell_A), 1);
        int solvertype = nrhs > 3 ? int(mxGetScalar(prhs[3])) : SLS_LU;
#pragma omp parallel for
        for (int i = 0; i < mxGetNumberOfElements(cell_A); i++)
        {
            const mxArray* mat = mxGetCell(cell_A, i);
            const mxArray* vec = mxGetCell(cell_b, i);
            mxArray* res = mxGetCell(plhs[0], i);
            bool suc = SparseLinearSolver::solve(solvertype, extractSpMat(mat), mxGetDoubles(vec), mxGetDoubles(res), mxGetN(vec));
            if (!suc) fprintf(stderr, "SparseLinearSolver failed to fullsolve\n");
            if (nlhs > 1) mxGetLogicals(plhs[1])[i] = suc;
        }
        return;
    }

    const pointer_t* ptr_handle = (pointer_t*)mxGetData(prhs[1]);
    std::vector<SparseLinearSolver*>* solver = (std::vector<SparseLinearSolver*>*)(*ptr_handle);
    for(int i = 0; i < solver->size(); i++)
        if (!solver->at(i)->isMemOk()) mexErrMsgTxt("bad memory! is the handle deleted already?");

    switch (mode) {
    case LST_CLOSE_HANDLE:
        {
        // deep free memory
        for (int i = 0; i < solver->size(); i++)
        {
            delete solver->at(i);
            solver->at(i) = 0;
        }
        delete solver;
        }
        break;
    case LST_NUMFACT:
    case LST_NUMFACT_AND_SOLVE:
    case LST_SOLVE:
        if (mode != LST_SOLVE) {
            // factorization only: SPLsolver(mode, handle, A, nonzeros); mode = 2 or 4
            if (nrhs < 4 + (mode == LST_NUMFACT_AND_SOLVE))
                mexErrMsgTxt("Invalid input: not enough input, missing new coefficient matrix (nonzero structure matrix + nonzero array) for factorization");

            //check input var
            const mxArray* cell_A = prhs[2];
            const mxArray* cell_nonzeros = prhs[3];
            if (!mxIsCell(cell_A) || !mxIsCell(cell_nonzeros))
                mexErrMsgTxt("Invalid input: A and nonzeros must be cells!");
            if (mxGetNumberOfElements(cell_A) != mxGetNumberOfElements(cell_nonzeros))
                mexErrMsgTxt("Invalid input: size(A) must equal to size(nonzeros)!");
            for (int i = 0; i < mxGetNumberOfElements(cell_A); i++)
            {
                const mxArray* mat_nonzeros = mxGetCell(cell_nonzeros, i);
                const mxArray* mat = mxGetCell(cell_A, i);
                if(!mxIsDouble(mat_nonzeros)) mexErrMsgTxt("Invalid input: nonzero vector should be in double precision.");
                if (mxGetNumberOfElements(mat_nonzeros) < nnz(mat)) mexErrMsgTxt("Invalid input: size of nonzero vector does not match the sparse matrix.");
            }
            
            //factorize
#pragma omp parallel for
            for (int i = 0; i < mxGetNumberOfElements(cell_A); i++)
            {
                const mxArray* mat_nonzeros = mxGetCell(cell_nonzeros, i);
                const mxArray* mat = mxGetCell(cell_A, i);
                bool suc = solver->at(i)->factorize(extractSpMat(mat, mxGetPr(mat_nonzeros)));
                if (!suc) fprintf(stderr, "SparseLinearSolver failed to factorize\n");
            }
        }

        if (mode == LST_NUMFACT) {
            plhs[0] = mxCreateLogicalMatrix(solver->size(), 1);
            for(int i = 0; i < solver->size(); i++)
                mxGetLogicals(plhs[0])[i] = mxCreateLogicalScalar(solver->at(i)->success());
        }
        else {
            // solve only: SPLsolver(3, handle, b);
            // numeric factorization and solve:  SPLsolver(4, handle, A, nonzeros, b)
            int iargb = 2 + (mode == LST_NUMFACT_AND_SOLVE) * 2;
            if (nrhs < iargb + 1)
                mexErrMsgTxt("Invalid input: not enough input, missing right hand side vector for solve");

            //check input var
            const mxArray* cell_b = prhs[iargb];
            if (!mxIsCell(cell_b))
                mexErrMsgTxt("Invalid input: b and nonzeros must be cells!");
            if (solver->size() != mxGetNumberOfElements(cell_b))
                mexErrMsgTxt("Invalid input: size(A) must equal to size(b)!");
            for (int i = 0; i < solver->size(); i++)
            {
                const mxArray* vec = mxGetCell(cell_b, i);
                if (!mxIsDouble(vec)) mexErrMsgTxt("Invalid input: vector b should be in double precision.");
                const size_t m = solver->at(i)->rows();
                if (mxGetM(vec) != m) mexErrMsgTxt("Invalid input: the size of vector b does not match the coefficient matrix.");
            }
            plhs[0] = mxDuplicateArray(cell_b);
            if (nlhs > 1) plhs[1] = mxCreateLogicalMatrix(solver->size(), 1);

#pragma omp parallel for
            for (int i = 0; i < solver->size(); i++)
            {
                const mxArray* vec = mxGetCell(cell_b, i);
                mxArray* res = mxGetCell(plhs[0], i);
                const size_t n = mxGetN(vec);
                bool suc = solver->at(i)->success() && solver->at(i)->solve(mxGetDoubles(vec), mxGetDoubles(res), n);
                if (!suc) fprintf(stderr, "SparseLinearSolver failed to solve\n");
                if (nlhs > 1) mxGetLogicals(plhs[1])[i] = suc;
            }          
        }
        break;

    default:
        mexErrMsgTxt("unrecognized mode");
    }
}
