/*!
 *  @file solver.hpp
 *  @brief These classes allow to solve a linear system.
 */

#ifndef SOLVER_HPP_
#define SOLVER_HPP_

#include "core/TypeDefinition.hpp"

//!
/*!
 * @class Solver
 * This is a base abstract class that implements a linear solver for the system Ax=b.
 * It receives a Eigen Sparse matrix A and a Eigen Vector b.
 * Each derived class employs a different solver.
 */
class Solver
{
public:

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
    Solver(const SpMat & A, const Vector & b):
        M_A(A), M_b(b) {}

    //! Get the matrix A
    /*!
     * @return the matrix A
     */
    const SpMat & getA() const
        { return M_A; }

    //! Get the RHS
    /*!
     * @return the vector b
     */
    const Vector & getb() const
        { return M_b; }

    //! Get the solution
    /*!
     * @return the vector containing the solution
     * @pre call solve()
     */
    const Vector & getSolution() const
        { return M_x; }

    //! Get the system size
    /*!
     * @return the system size, i.e. the number of rows of A
     */
    UInt size() const
        { return M_A.rows(); }

    //! Get the non zero values
    /*!
     * @return the number or non zero values
     */
    UInt nonZeros() const
        { return M_A.nonZeros(); }

    //! Solve the linear system
    /*!
     * Solve the linear system Ax=b.
     * The method employed to solve the system depends on the derived class of Solver that calls this method.
     */
    virtual void solve() = 0;

    //! Destructor
    virtual ~Solver() {}

protected:

    //! Sparse Matrix
    const SpMat & M_A;
    //! RHS
    const Vector & M_b;
    //! Solution vector
    Vector M_x;

};

class EigenCholesky : public Solver
{
public:

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
    EigenCholesky(const SpMat & A, const Vector & b):
        Solver(A,b) {}

    //! Solve the linear system
    /*!
     * Solve the linear system Ax=b by means of a Cholesky factorization.
     */
    virtual void solve();

    //! Destructor
    virtual ~EigenCholesky() {}

};

class EigenLU : public Solver
{
public:

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
    EigenLU(const SpMat & A, const Vector & b):
        Solver(A,b) {}

    //! Solve the linear system
    /*!
     * Solve the linear system Ax=b by means of a Cholesky factorization.
     */
    virtual void solve();

    //! Destructor
    virtual ~EigenLU() {}

};

#endif /* SOLVER_HPP_ */
