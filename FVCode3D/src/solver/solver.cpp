/*!
 *  @file solver.cpp
 *  @brief These classes allow to solve a linear system (definitions).
 */

#include "solver/solver.hpp"

void EigenCholesky::solve()
{
    Eigen::SimplicialCholesky<SpMat> chol(*M_A);
    M_x = chol.solve(*M_b);
}

void EigenLU::solve()
{
    Eigen::SparseLU<SpMat> lu;
    lu.analyzePattern(*M_A);
    lu.factorize(*M_A);
    M_x = lu.solve(*M_b);
} // solve

void EigenUmfPack::solve()
{
    Eigen::UmfPackLU<SpMat> lu( *M_A );
    M_x = lu.solve( *M_b );
} // solve

void EigenCG::solve()
{
	Eigen::ConjugateGradient<SpMat> cg;
	
	cg.setMaxIterations(M_maxIter);
	cg.setTolerance(M_tol);
	
	cg.compute(*M_A);
	M_x = cg.solve(*M_b);
	
	M_iter = cg.iterations();
	M_res = cg.error();
}

void EigenBiCGSTAB::solve()
{
	Eigen::BiCGSTAB<SpMat> bicgstab;
	
	bicgstab.setMaxIterations(M_maxIter);
	bicgstab.setTolerance(M_tol);
	
	bicgstab.compute(*M_A);
	M_x = bicgstab.solve(*M_b);
	
	M_iter = bicgstab.iterations();
	M_res = bicgstab.error();	
}
