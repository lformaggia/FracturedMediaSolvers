/*!
 *	@file solver.cpp
 *	@brief These classes allow to solve a linear system (definitions).
 */

#include "solver/solver.hpp"

void EigenCholesky::solve()
{
	Eigen::SimplicialCholesky<SpMat> chol(M_A);
	M_x = chol.solve(M_b);
}
