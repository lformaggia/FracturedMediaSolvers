/*!
 *  @file solver.cpp
 *  @brief These classes allow to solve a linear system (definitions).
 */

#include "solver/solver.hpp"

void EigenCholesky::solve()
{
    Eigen::SimplicialCholesky<SpMat> chol(M_A);
    M_x = chol.solve(M_b);
}

void EigenLU::solve()
{
    Eigen::SparseLU<SpMat> lu;
    lu.analyzePattern(M_A);
    lu.factorize(M_A);
    M_x = lu.solve(M_b);
} // solve
