/*!
 * @file solver.hpp
 * @brief These classes allow to solve a linear system.
 */

#ifndef SOLVER_HPP_
#define SOLVER_HPP_

#include <FVCode3D/core/TypeDefinition.hpp>

namespace FVCode3D
{

class Solver;

//! Typedef for SolverPtr_Type
/*!
 * @typedef SolverPtr_Type
 * This type definition permits to handle a SolverPtr_Type as a SolverPtr_Type.
 */
typedef std::shared_ptr<Solver> SolverPtr_Type;

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
    Solver( const UInt nbDofs = 0 ) :
        M_A( nbDofs, nbDofs ),
        M_b( Vector::Zero( nbDofs ) ),
        M_x( Vector::Zero( nbDofs ) )
    {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
    Solver(const SpMat & A, const Vector & b):
        M_A(A), M_b(b) {}

    //! Set the matrix A
    /*!
     * @param the matrix A
     */
    void setA(const SpMat & A) { M_A = A; }

    //! Set the RHS
    /*!
     * @param the vector b
     */
    void setb(const Vector & b) { M_b = b; }

    //! Resize the matrix and vector
    /*!
     * @param nbDofs number of the dofs
     */
    void setDofs(const UInt nbDofs)
    {
        M_A.resize(nbDofs, nbDofs);
        M_b = Vector::Zero( nbDofs );
        M_x = Vector::Zero( nbDofs );
    }

    //! Get the matrix A
    /*!
     * @return the matrix A
     */
    const SpMat & getA() const { return M_A; }

    //! Get the matrix A
    /*!
     * @return the matrix A
     */
    SpMat & getA() { return M_A; }

    //! Get the RHS
    /*!
     * @return the vector b
     */
    const Vector & getb() const { return M_b; }

    //! Get the RHS
    /*!
     * @return the vector b
     */
    Vector & getb() { return M_b; }

    //! Get the solution
    /*!
     * @return the vector containing the solution
     * @pre call solve()
     */
    const Vector & getSolution() const { return M_x; }

    //! Get the solution
    /*!
     * @return the vector containing the solution
     * @pre call solve()
     */
    Vector & getSolution() { return M_x; }

    //! Get the system size
    /*!
     * @return the system size, i.e. the number of rows of A
     */
    UInt size() const { return M_A.rows(); }

    //! Get the non zero values
    /*!
     * @return the number or non zero values
     */
    UInt nonZeros() const { return M_A.nonZeros(); }

    //! Solve the linear system
    /*!
     * Solve the linear system Ax=b.
     * The method employed to solve the system depends on the derived class of Solver that calls this method.
     */
    virtual void solve() = 0;

    //! Destructor
    virtual ~Solver() = default;

protected:

    //! Sparse Matrix
    SpMat M_A;

    //! RHS
    Vector M_b;

    //! Solution vector
    Vector M_x;
}; // class Solver

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
    EigenCholesky( const UInt nbDofs = 0 ): Solver( nbDofs ) {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
    EigenCholesky(const SpMat & A, const Vector & b):
        Solver(A, b) {}

    //! Solve the linear system
    /*!
     * Solve the linear system Ax=b by means of a Cholesky factorization.
     */
    virtual void solve();

    //! Destructor
    virtual ~EigenCholesky() = default;
}; // class EigenCholesky

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
    EigenLU( const UInt nbDofs = 0 ): Solver( nbDofs ) {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
    EigenLU(const SpMat & A, const Vector & b):
        Solver(A, b) {}

    //! Solve the linear system
    /*!
     * Solve the linear system Ax=b by means of a LU factorization.
     */
    virtual void solve();

    //! Destructor
    virtual ~EigenLU() = default;
}; // class EigenLU

#ifdef FVCODE3D_HAS_UMFPACK
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
    EigenUmfPack( const UInt nbDofs = 0 ): Solver( nbDofs ) {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
    EigenUmfPack(const SpMat & A, const Vector & b):
        Solver(A, b) {}

    //! Solve the linear system
    /*!
     * Solve the linear system Ax=b by means of a UmfPack factorization.
     */
    virtual void solve();

    //! Destructor
    virtual ~EigenUmfPack() = default;
};
#endif // FVCODE3D_HAS_UMFPACK


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
    IterativeSolver( const UInt nbDofs = 0 ):
        Solver( nbDofs ),
        M_maxIter( S_referenceMaxIter ), M_iter(0), M_tol( S_referenceTol ), M_res(0), CIndex(0) {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
    IterativeSolver( const SpMat & A, const Vector & b,
                     const UInt maxIter = S_referenceTol, const Real tol = S_referenceMaxIter ):
        Solver( A, b), M_maxIter(maxIter), M_iter(0), M_tol(tol), M_res(0), CIndex(0) {}

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

    //! Get the normalized residual error
    /*!
     * @return the normalized residual error || b - A * x || / || b ||
     */
    Real getResidual() const { return M_res; }
    //@}

    //! @name Set Methods
    //@{

    //! Set maximum number of iterations
    /*!
     * @param maximum number of iterations
     */
    void setMaxIter(const UInt maxIter) { M_maxIter = M_iter = maxIter; };

    //! Set relative tolerance for the residual error
    /*!
     * @param relative tolerance for the residual error
     */
    void setTolerance(const Real tol) { M_tol = M_res = tol; };
    //@}
    
    //! Print out the computation details
    /*!
     * Print out the cimputation details
     */
    virtual void print()
    {
		std::cout << std::endl;
		std::cout << "\t# iterations: " << M_iter << std::endl;
		std::cout << "\tResidual: " << M_res << std::endl<<std::endl;
		switch(CIndex)
		{
			case 0 :
			std::cout<<"Ok, the method converge within the tolerance."<<std::endl;
			break;
			case 1 :
			std::cout<<"The method does not converge."<<std::endl;
		}
	};

    //! Solve the linear system
    /*!
     * Solve the linear system Ax=b.
     * The method employed to solve the system depends on the derived class of Solver that calls this method.
     */
    virtual void solve() = 0;

    //! Destructor
    virtual ~IterativeSolver() = default;

protected:

    //! Maximum number of iterations
    UInt M_maxIter;
    //! Number of iterations
    UInt M_iter;
    //! Tolerance
    Real M_tol;
    //! Residual error
    Real M_res;
    //! Computation index
    UInt CIndex;

    static constexpr Real S_referenceTol = 1e-6;
    static constexpr UInt S_referenceMaxIter = 100;
}; // class IterativeSolver

//! Class imlCG
/*!
 * @class imlCG
 * This class implements a linear solver for the system Ax=b.
 * It uses the conjugate gradient method on a spd matrix.
 */
class imlCG : public IterativeSolver
{
public:

    //! Empty constructor
    imlCG( const UInt nbDofs = 0 ): IterativeSolver( nbDofs ) {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
    imlCG( const SpMat & A, const Vector & b,
             const UInt maxIter = IterativeSolver::S_referenceTol,
             const Real tol = IterativeSolver::S_referenceMaxIter ):
        IterativeSolver(A, b,maxIter,tol) {}
        
    //! Destructor
    virtual ~imlCG() = default;

    //! Solve the linear system
    /*!
     * Solve the linear system Ax=b by means of the conjugate gradient method.
     */
    virtual void solve();
    
}; // class EigenCG


//! Class imlBiCGSTAB
/*!
 * @class imlBiCGSTAB
 * This class implements a linear solver for the system Ax=b.
 * It uses the stabilized bi-conjugate gradient method on a square matrix.
 */
class imlBiCGSTAB : public IterativeSolver
{
public:

    //! Empty constructor
    imlBiCGSTAB( const UInt nbDofs = 0 ): IterativeSolver( nbDofs ), restart(Default_restart) {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
    imlBiCGSTAB( const SpMat & A, const Vector & b,
				const UInt maxIter = IterativeSolver::S_referenceTol,
				const Real tol = IterativeSolver::S_referenceMaxIter, const bool rest = Default_restart ):
        IterativeSolver(A, b,maxIter,tol), restart(rest) {}
        
    //! Destructor
    virtual ~imlBiCGSTAB() = default;

    //! Set the restart
    /*!
     * @param The restart
     */
    void setRestart(bool rest){ restart = rest; };
    
    //! Print out the computation details
    /*!
     * Print out the cimputation details
     */
    void print()
    {
		IterativeSolver::print();
		switch(CIndex)
		{
			case 2 :
			std::cout<<"A 1st type breakdown is occurred."<<std::endl;
			break;
			case 3 :
			std::cout<<"A 2st type breakdown is occurred."<<std::endl;
		}
	};

    //! Solve the linear system
    /*!
     * Solve the linear system Ax=b by means of the stabilized bi-conjugate gradient method.
     */
    virtual void solve();
    
private:
	//! A bool to know if using BiCGSTAB with restart or not
	bool restart;
	//! The default restart value
	static constexpr bool Default_restart = false;
	
}; // class BiCGSTAB


//! Class imlGMRES
/*!
 * @class imlGMRES
 * This class implements a linear solver for the system Ax=b.
 * It uses the generalized minimum residual method with restart on a square matrix.
 */
class imlGMRES : public IterativeSolver
{
public:

    //! Empty constructor
    imlGMRES( const UInt nbDofs = 0 ): IterativeSolver( nbDofs ), m(Default_m) {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
    imlGMRES( const SpMat & A, const Vector & b,
				const UInt maxIter = IterativeSolver::S_referenceTol,
				const Real tol = IterativeSolver::S_referenceMaxIter, const UInt rest = Default_m ):
        IterativeSolver(A, b,maxIter,tol), m(rest) {}
  
    //! Destructor
    virtual ~imlGMRES() = default;
          
    //! Set the restart
    /*!
     * @param The restart
     */
    void set_m(UInt rest_val){ m = rest_val; };
    
    //! Solve the linear system
    /*!
     * Solve the linear system Ax=b by means of the GMRES method.
     */
    virtual void solve();

private:
	//! The restart level (how many iterations nedded to perform a restart)
	UInt m;
	//! The default restart value
	static constexpr UInt Default_m = 40;
	
}; // class GMRES


#ifdef FVCODE3D_HAS_SAMG
//! Class SamgSolver
/*!
 * @class SamgSolver
 * This class implements a linear solver for the system Ax=b.
 * It uses algebraic multi-grid methods (SAMG library).
 */
class SamgSolver : public IterativeSolver
{
public:

    //! Empty constructor
    SamgSolver( const UInt nbDofs = 0 ): IterativeSolver( nbDofs ) {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
    SamgSolver( const SpMat & A, const Vector & b,
                const UInt maxIter = IterativeSolver::S_referenceTol,
                const Real tol = IterativeSolver::S_referenceMaxIter ):
        IterativeSolver(A, b,maxIter,tol) {}

    //! Solve the linear system
    virtual void solve();

    //! Destructor
    virtual ~SamgSolver() = default;

protected:

    class SamgParameters;

    //! Set the parameters for the specific method
    /*!
     * @param SP SamgParameters instance
     */
    virtual void setSamgParameters(SamgParameters & SP) = 0;

    //! Class SamgParameters
    /*!
     * @class SamgParameters
     * This class contains the parameters used to the SAMG solver.
     */
    class SamgParameters
    {
    public:

        //! Default constructor
        /*!
         * Initialize the parameters
         */
        SamgParameters();

        //! No copy constructor
        SamgParameters(const SamgParameters &) = delete;

        //! Destructor
        ~SamgParameters();

    public:

        //! Primary parameters
        //@{
        //! Dimension of (dummy) vector iu
        int    ndiu;
        //! Dimension of (dummy) vector ip
        int    ndip;
        //! Solution approach
        int    nsolve;
        //! First approximation of the solution
        int    ifirst;
        //! Required residual reduction
        double eps;
        //! Cycling and acceleration strategy
        int    ncyc;
        //! Estimated dimensioning: operator completity
        double a_cmplx;
        //! Estimated dimensioning: grid complexity
        double g_cmplx;
        //! Estimated dimensioning: average row length of interpolation
        double w_avrge;
        //! Estimated dimensioning: point complexity
        double p_cmplx;
        //! Input matrix checking
        double chktol;
        //! Printed output during the setup phase
        int    idump;
        //! Printed output during the solution phase
        int    iout;
        //! Select coarsening strategy
        int    n_default;
        //! General control switch
        int    iswtch;
        //@}

        //! AMG declarations
        //@{
            //! Input
            //@{
        //! Number of points
        int npnt;
        //! Number of unknown
        int nsys;
        //! Matrix type
        int matrix;
        //! Number of variables
        int nnu;
        //! Number of entries
        int nna;
        //! Map variable to unknown
        int * iu;
        //! Map variable to point
        int * ip;
        //! Scale the output solution
        int * iscale;
            //@}

            //! Output
            //@{
        //! Code indicating warnings(<0) or errors(>0)
        int ierr;
        //! Code indicating warnings(<0) or errors(>0) at samg_leave call
        int ierrl;
        //! Total number of cycles
        int ncyc_done;
        //! Residual of final approximation
        double res_out;
        //! Residual of first guess
        double res_in;
            //@}
        //@}
    }; // class SamgParameters
}; // class SamgSolver

//! Class SamgSym
/*!
 * @class SamgSym
 * This class implements a linear solver for the system Ax=b.
 * It uses algebraic multi-grid methods (SAMG library).
 * It is used in the case of symmetric matrix A.
 */
class SamgSym : public SamgSolver
{
public:

    //! Empty constructor
    SamgSym( const UInt nbDofs = 0 ): SamgSolver( nbDofs ) {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
    SamgSym( const SpMat & A, const Vector & b,
             const UInt maxIter = IterativeSolver::S_referenceTol,
             const Real tol = IterativeSolver::S_referenceMaxIter ):
        SamgSolver(A, b,maxIter,tol) {}

    //! Destructor
    virtual ~SamgSym() = default;

protected:

    //! Set the parameters for the specific method
    /*!
     * @param SP SamgParameters instance
     */
    virtual void setSamgParameters(SamgSolver::SamgParameters & SP);

}; // class SamgSym

//! Class SamgNotSym
/*!
 * @class SamgNotSym
 * This class implements a linear solver for the system Ax=b.
 * It uses algebraic multi-grid methods (SAMG library).
 * It is used in the case of not symmetric matrix A.
 */
class SamgNotSym : public SamgSolver
{
public:

    //! Empty constructor
    SamgNotSym( const UInt nbDofs = 0 ): SamgSolver( nbDofs ) {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
    SamgNotSym( const SpMat & A, const Vector & b,
                const UInt maxIter = IterativeSolver::S_referenceTol,
                const Real tol = IterativeSolver::S_referenceMaxIter ):
        SamgSolver(A, b,maxIter,tol) {}

    //! Destructor
    virtual ~SamgNotSym() = default;

protected:

    //! Set the parameters for the specific method
    /*!
     * @param SP SamgParameters instance
     */
    virtual void setSamgParameters(SamgSolver::SamgParameters & SP);

}; // class SamgNotSym
#endif // FVCODE3D_HAS_SAMG

} // namespace FVCode3D
#endif /* SOLVER_HPP_ */
