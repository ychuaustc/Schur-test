#ifndef SPLSOLVER_H
#define SPLSOLVER_H
#include <vector>
#include <variant>

#include <Eigen/Dense>
#include <Eigen/SparseCore>

#include<Eigen/SparseCholesky>
#include<Eigen/SparseLU>


using index_t = int64_t;
//using index_t = int;
using SpMat = Eigen::SparseMatrix<double, 0, index_t>;

enum SPLINEAR_SOLVER_DIRECT {
    SLS_LLT = 0, SLS_LDLT, SLS_LU
};

struct SparseLinearSolver
{
private:
    const int SLS_INTEGRITY_CHECK = 0x12345678;
    int memIntFlag = SLS_INTEGRITY_CHECK;

public:
    using solver_t = std::variant<Eigen::SimplicialLLT<SpMat>, Eigen::SimplicialLDLT<SpMat>, Eigen::SparseLU<SpMat>>;
    solver_t imp;

    SparseLinearSolver(int s = 0) { emplaceSolver(imp, s); }

    static void emplaceSolver(solver_t& solver, int s)
    {
        switch (s) {
        case SLS_LLT: solver.emplace<SLS_LLT>(); break;
        case SLS_LDLT: solver.emplace<SLS_LDLT>(); break;
        case SLS_LU: solver.emplace<SLS_LU>(); break;
        default:
            fprintf(stderr, "unexpected linear solver %d!", s);
        }
    }

    template<class mat>
    static bool solve(int s, const mat& M, const double* b, double* x, size_t n) {
        const size_t m = M.rows();
        solver_t solver;
        emplaceSolver(solver, s);
        return std::visit([&](auto& s) {
            s.compute(M);
            if (s.info() != Eigen::Success) return false;
            Eigen::MatrixXd tmp;
            tmp.resize(m, n);
            s._solve_impl(Eigen::Map<const Eigen::MatrixXd>(b, m, n), tmp); //why Eigen::Map<Eigen::MatrixXd> can't convert to Eigen::MatrixXd& ?
            Eigen::Map<Eigen::MatrixXd>(x, m, n) = tmp;
            return s.info() == Eigen::Success;
            }, solver);
    }

    inline bool solve(const double* b, double* x, size_t n) {
        return std::visit([&](auto& solver) {
            const size_t m = solver.cols();
            Eigen::MatrixXd tmp;
            tmp.resize(m, n);
            solver._solve_impl(Eigen::Map<const Eigen::MatrixXd>(b, m, n), tmp);
            Eigen::Map<Eigen::MatrixXd>(x, m, n) = tmp;
            return solver.info() == Eigen::Success;
            }, imp);
    }

    // use template matrix type instead of SpMat, otherwise a Map<SpMat> will be implicitly converted to a SpMat, which will be destructed later and causing problem
    // not sure if the solver will access the nonzeros during symbolic factorization
    // inline void analyzePattern(const SpMat& M) 
    template<class mat>
    inline void analyzePattern(const mat& M) { std::visit([&](auto& solver) { solver.analyzePattern(M); }, imp); }

    template<class mat>
    inline void compute(const mat& M) { std::visit([&](auto& solver) { solver.compute(M); }, imp); }

    template<class mat>
    inline bool factorize(const mat& M) {
        return std::visit([&](auto& solver) {
            solver.factorize(M);
            return solver.info() == Eigen::Success;
            }, imp);
    }

    inline size_t rows() { return std::visit([&](const auto& solver) { return solver.rows(); }, imp); }

    inline bool success() { return std::visit([&](const auto& solver) { return solver.info() == Eigen::Success; }, imp); }

    ~SparseLinearSolver() { memIntFlag = 0; /*cancel the integrity check*/ }
    bool isMemOk() const { return memIntFlag == SLS_INTEGRITY_CHECK; }
};

#endif
