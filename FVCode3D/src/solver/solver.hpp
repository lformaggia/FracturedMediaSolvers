/*!
 *  @file solver.hpp
 *  @brief These classes allow to solve a linear system.
 */

#ifndef SOLVER_HPP_
#define SOLVER_HPP_

#include "core/TypeDefinition.hpp"

//! Class Solver
/*!
 * @class Solver
 * This is a base abstract class that implements a linear solver for the system Ax=b.
 * It receives a Eigen Sparse matrix A and a Eigen Vector b.
 * Each derived class employs a different solver.
 */
class Solver
{
public:

    //! Empty constructor
    Solver():M_A(nullptr), M_b(nullptr) {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
    Solver(const SpMat & A, const Vector & b):
        M_A(&A), M_b(&b) {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
    Solver(const SpMat * A, const Vector * b):
        M_A(A), M_b(b) {}

    //! Set the matrix A
    /*!
     * @param the matrix A
     */
    void setA(const SpMat & A)
    	{ M_A = &A; }

    //! Set the RHS
    /*!
     * @param the vector b
     */
    void setb(const Vector & b)
        { M_b = &b; }

    //! Get the matrix A
    /*!
     * @return the matrix A
     */
    const SpMat & getA() const
        { return *M_A; }

    //! Get the RHS
    /*!
     * @return the vector b
     */
    const Vector & getb() const
        { return *M_b; }

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
        { return M_A->rows(); }

    //! Get the non zero values
    /*!
     * @return the number or non zero values
     */
    UInt nonZeros() const
        { return M_A->nonZeros(); }

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
    const SpMat * M_A;
    //! RHS
    const Vector * M_b;
    //! Solution vector
    Vector M_x;

};

//! Class EigenCholesky
/*!
 * @class EigenCholesky
 * This class implements a linear solver for the system Ax=b.
 * It uses the Cholesky factorization on a SPD matrix.
 */
class EigenCholesky : public Solver
{
public:

    //! Empty constructor
	EigenCholesky(): Solver() {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
    EigenCholesky(const SpMat & A, const Vector & b):
        Solver(&A,&b) {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
    EigenCholesky(const SpMat * A, const Vector * b):
        Solver(A,b) {}

    //! Solve the linear system
    /*!
     * Solve the linear system Ax=b by means of a Cholesky factorization.
     */
    virtual void solve();

    //! Destructor
    virtual ~EigenCholesky() {}

};

//! Class EigenLU
/*!
 * @class EigenLU
 * This class implements a linear solver for the system Ax=b.
 * It uses the LU factorization on a square matrix.
 */
class EigenLU : public Solver
{
public:

    //! Empty constructor
	EigenLU(): Solver() {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
    EigenLU(const SpMat & A, const Vector & b):
        Solver(&A,&b) {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
    EigenLU(const SpMat * A, const Vector * b):
        Solver(A,b) {}

    //! Solve the linear system
    /*!
     * Solve the linear system Ax=b by means of a LU factorization.
     */
    virtual void solve();

    //! Destructor
    virtual ~EigenLU() {}

};

//! Class EigenUmfPack
/*!
 * @class EigenUmfPack
 * This class implements a linear solver for the system Ax=b.
 * It uses the LU factorization on a square matrix (from UmfPack).
 */
class EigenUmfPack : public Solver
{
public:

	//! Empty constructor
	EigenUmfPack(): Solver() {}

	//! Constructor
	/*!
	 * @param A Eigen sparse matrix
	 * @param b RHS, it is Eigen vector
	 */
	EigenUmfPack(const SpMat & A, const Vector & b):
		Solver(&A,&b) {}

	//! Constructor
	/*!
	 * @param A Eigen sparse matrix
	 * @param b RHS, it is Eigen vector
	 */
	EigenUmfPack(const SpMat * A, const Vector * b):
		Solver(A,b) {}

	//! Solve the linear system
	/*!
	 * Solve the linear system Ax=b by means of a UmfPack factorization.
	 */
	virtual void solve();

	//! Destructor
	virtual ~EigenUmfPack() {}
};

//! Class IterativeSolver
/*!
 * @class IterativeSolver
 * This is an abstract class that implements a linear solver for the system Ax=b solved by means of an iterative method.
 * It receives a Eigen Sparse matrix A and a Eigen Vector b.
 * Each derived class employs a different solver.
 */
class IterativeSolver : public Solver
{
public:

	//! Empty constructor
	IterativeSolver(): Solver(), M_maxIter(100), M_iter(0), M_tol(1e-6), M_res(0) {}

	//! Constructor
	/*!
	 * @param A Eigen sparse matrix
	 * @param b RHS, it is Eigen vector
	 */
	IterativeSolver(const SpMat & A, const Vector & b):
		Solver(&A,&b), M_maxIter(100), M_iter(0), M_tol(1e-6), M_res(0) {}

	//! Constructor
	/*!
	 * @param A Eigen sparse matrix
	 * @param b RHS, it is Eigen vector
	 */
	IterativeSolver(const SpMat * A, const Vector * b):
		Solver(A,b), M_maxIter(100), M_iter(0), M_tol(1e-6), M_res(0) {}

	//! Constructor
	/*!
	 * @param A Eigen sparse matrix
	 * @param b RHS, it is Eigen vector
	 */
	IterativeSolver(const SpMat & A, const Vector & b, const UInt maxIter, const Real tol):
		Solver(&A,&b), M_maxIter(maxIter), M_iter(0), M_tol(tol), M_res(0) {}

	//! Constructor
	/*!
	 * @param A Eigen sparse matrix
	 * @param b RHS, it is Eigen vector
	 */
	IterativeSolver(const SpMat * A, const Vector * b, const UInt maxIter, const Real tol):
		Solver(A,b), M_maxIter(maxIter), M_iter(0), M_tol(tol), M_res(0) {}

	//! @name Get Methods
	//@{

	//! Get maximum number of iterations
	/*!
	 * @return maximum number of iterations
	 */
	UInt getMaxIter() const { return M_maxIter; }

	//! Get number of iterations
	/*!
	 * @return number of iterations employed by the solver
	 */
	UInt getIter() const { return M_iter; }

	//! Get relative tolerance for the residual error
	/*!
	 * @return relative tolerance for the residual error
	 */
	Real getTolerance() const { return M_tol; }

	//! Get the residual error
	/*!
	 * @return residual error ( b - A * x )
	 */
	Real getResidual() const { return M_res; }

	//@}

	//! @name Set Methods
	//@{

	//! Set maximum number of iterations
	/*!
	 * @param maximum number of iterations
	 */
	void setMaxIter(const UInt maxIter) { M_maxIter = maxIter; };

	//! Set relative tolerance for the residual error
	/*!
	 * @param relative tolerance for the residual error
	 */
	void setTolerance(const Real tol) { M_tol = tol; };

	//@}

	//! Solve the linear system
	/*!
	 * Solve the linear system Ax=b.
	 * The method employed to solve the system depends on the derived class of Solver that calls this method.
	 */
	virtual void solve() = 0;

	//! Destructor
	virtual ~IterativeSolver() {}

protected:

	//! Maximum number of iterations
	UInt M_maxIter;
	//! Number of iterations
	UInt M_iter;
	//! Tolerance
	Real M_tol;
	//! Residual error
	Real M_res;

};

//! Class EigenCG
/*!
 * @class EigenCG
 * This class implements a linear solver for the system Ax=b.
 * It uses the conjugate gradient method on a spd matrix.
 */
class EigenCG : public IterativeSolver
{
public:

    //! Empty constructor
	EigenCG(): IterativeSolver() {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
	EigenCG(const SpMat & A, const Vector & b):
		IterativeSolver(&A,&b) {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
	EigenCG(const SpMat * A, const Vector * b):
		IterativeSolver(A,b) {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
	EigenCG(const SpMat & A, const Vector & b, const UInt maxIter, const Real tol):
		IterativeSolver(&A,&b,maxIter,tol) {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
	EigenCG(const SpMat * A, const Vector * b, const UInt maxIter, const Real tol):
		IterativeSolver(A,b,maxIter,tol) {}

    //! Solve the linear system
    /*!
     * Solve the linear system Ax=b by means of the conjugate gradient method.
     */
    virtual void solve();

    //! Destructor
    virtual ~EigenCG() {}

};

//! Class EigenBiCGSTAB
/*!
 * @class EigenBiCGSTAB
 * This class implements a linear solver for the system Ax=b.
 * It uses the stabilized bi-conjugate gradient method on a square matrix.
 */
class EigenBiCGSTAB : public IterativeSolver
{
public:

    //! Empty constructor
	EigenBiCGSTAB(): IterativeSolver() {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
	EigenBiCGSTAB(const SpMat & A, const Vector & b):
		IterativeSolver(&A,&b) {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
	EigenBiCGSTAB(const SpMat * A, const Vector * b):
		IterativeSolver(A,b) {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
	EigenBiCGSTAB(const SpMat & A, const Vector & b, const UInt maxIter, const Real tol):
		IterativeSolver(&A,&b,maxIter,tol) {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
	EigenBiCGSTAB(const SpMat * A, const Vector * b, const UInt maxIter, const Real tol):
		IterativeSolver(A,b,maxIter,tol) {}	

    //! Solve the linear system
    /*!
     * Solve the linear system Ax=b by means of the conjugate gradient method.
     */
    virtual void solve();

    //! Destructor
    virtual ~EigenBiCGSTAB() {}

};

#endif /* SOLVER_HPP_ */
