/*!
 *  @file solver.cpp
 *  @brief These classes allow to solve a linear system (definitions).
 */

#include <FVCode3D/solver/Solver.hpp>
#ifdef FVCODE3D_HAS_UMFPACK
#include <Eigen/UmfPackSupport>
#endif

namespace FVCode3D
{

void EigenCholesky::solve()
{
    Eigen::SimplicialCholesky<SpMat> chol(M_A);
    M_x = chol.solve(M_b);
} // EigenCholesky::solve

void EigenLU::solve()
{
    Eigen::SparseLU<SpMat> lu;
    lu.analyzePattern(M_A);
    lu.factorize(M_A);
    M_x = lu.solve(M_b);
} // EigenLU::solve

#ifdef FVCODE3D_HAS_UMFPACK
void EigenUmfPack::solve()
{
    Eigen::UmfPackLU<SpMat> lu( M_A );
    M_x = lu.solve( M_b );
} // EigenUmfPack::solve
#endif // FVCODE3D_HAS_UMFPACK

Real IterativeSolver::S_referenceTol = 1e-6;
UInt IterativeSolver::S_referenceMaxIter = 100;

void EigenCG::solve()
{
    Eigen::ConjugateGradient<SpMat> cg;

    cg.setMaxIterations(M_maxIter);
    cg.setTolerance(M_tol);

    cg.compute(M_A);
    M_x = cg.solve(M_b);

    M_iter = cg.iterations();
    M_res = cg.error();
} // EigenCG::solve

void EigenBiCGSTAB::solve()
{
    Eigen::BiCGSTAB<SpMat> bicgstab;

    bicgstab.setMaxIterations(M_maxIter);
    bicgstab.setTolerance(M_tol);

    bicgstab.compute(M_A);
    M_x = bicgstab.solve(M_b);

    M_iter = bicgstab.iterations();
    M_res = bicgstab.error();
} // EigenBiCGSTAB::solve

} // namespace FVCode3D
