/*!
 * @file solver.hpp
 * @brief These classes allow to solve a linear system.
 */

#ifndef SOLVER_HPP_
#define SOLVER_HPP_

#include <string>

#include <FVCode3D/core/TypeDefinition.hpp>
#include <FVCode3D/preconditioner/preconditioner.hpp>
#include <FVCode3D/preconditioner/preconHandler.hpp>

namespace FVCode3D
{
	
class Solver;
	
//! Typedef for SolverPtr_Type
/*!
 * @typedef SolverPtr_Type
 * This type definition permits to handle a std::shared_ptr<Solver> as a SolverPtr_Type.
 */
typedef std::shared_ptr<Solver> SolverPtr_Type;

//! Class Solver
/*!
 * @class Solver
 * This is a base abstract class that implements a linear solver for the system Ax=b.
 * The derived abstract classes are DirectSolver and IterativeSolver, the first stores 
 * the monolithic matrix, the second stores the single block matrices.
 */
class Solver
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Empty constructor
    Solver(const UInt nbDofs = 0):
		M_b(Vector::Zero(nbDofs)),
		M_x(Vector::Zero(nbDofs)) {}

    //! Constructor
    /*!
     * @param b RHS, it is Eigen vector
     */
    Solver(const Vector & b):
        M_b(b) {}
    
    //! Destructor
    virtual ~Solver() = default;
    //@}

    //! @name Set Methods
    //@{
    //! Set the RHS
    /*!
     * @param the vector b
     */
    void setb(const Vector & b) { M_b = b; }

    //! Set the system dofs
    /*!
     * @param nbDofs number of the dofs
     */
    virtual void setDofs(const UInt nbDofs)
    {
        M_b = Vector::Zero( nbDofs );
        M_x = Vector::Zero( nbDofs );
    }
    //@}

    //! @name Get Methods
    //@{
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
    
    //! Size of the system
    /*!
     * @return The system size
     */   
    UInt getSize() const { return M_b.size(); }
    //@}

	//! @name Methods
    //@{   
    //! Solve the linear system
    /*!
     * Solve the linear system Ax=b.
     * The method employed to solve the system depends on the derived class of Solver that calls this method.
     */
    virtual void solve() = 0;
    //@}

protected:
    //! RHS
    Vector      		M_b;
    //! Solution vector
    Vector      		M_x;
}; // class Solver

//! Class DirectSolver
/*!
 * @class DirectSolver
 * It implements the base class for a direct solver like, for example, the LU solver or Cholesky solver.
 * The matrix stored is that of the monolithic system that has the type of an Eigen Sparse Matrix.
 * Each derived class employes a different direct solver.
 */
class DirectSolver : public Solver
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Empty constructor
    DirectSolver(const UInt nbDofs = 0): 
		M_A(nbDofs,nbDofs),
		Solver(nbDofs) {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
    DirectSolver( const SpMat & A, const Vector & b):
        M_A(A), Solver(b) {}
        
    //! Destructor
    virtual ~DirectSolver() = default;
	//@}

    //! @name Set Methods
    //@{
    //! Set the system matrix
    /*!
     * @param A the system matrix
     */
    void setA(const SpMat & A) { M_A = A; }

    //! Set the system dofs
    /*!
     * @param nbDofs number of the dofs
     */
    void setDofs(const UInt nbDofs)
    {
		M_A.resize(nbDofs, nbDofs);
		Solver::setDofs(nbDofs);
    }
    //@}

    //! @name Get Methods
    //@{  
    //! Get the system matrix
    /*!
     * @return the system matrix
     */
    const SpMat & getA() const { return M_A; }

    //! Get the system matrix
    /*!
     * @return the system matrix
     */
    SpMat & getA() { return M_A; }
    
    //! Get the number of non zero
    /*!
     * @return the number of non zero
     */
    UInt getNonZero() { return M_A.nonZeros(); }
    //@}

    //! @name Methods
    //@{
    //! Solve the linear system
    /*!
     * Solve the linear system Ax=b.
     * The method employed to solve the system depends on the derived class of Solver that calls this method.
     */
    virtual void solve() = 0;
    //@}
    
protected:
	//! The system matrix
	SpMat    M_A;
    
}; // class DirectSolver

//! Class EigenCholesky
/*!
 * @class EigenCholesky
 * This class implements a linear solver for the system Ax=b.
 * It uses the Cholesky factorization on a SPD matrix.
 */
class EigenCholesky : public DirectSolver
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Empty constructor
    EigenCholesky(const UInt nbDofs = 0): 
		DirectSolver(nbDofs) {}
    
    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
    EigenCholesky(const SpMat & A, const Vector & b):
        DirectSolver(A, b) {}
    
    //! Destructor
    ~EigenCholesky() = default;   
    //@}

    //! @name Methods
    //@{
    //! Solve the linear system
    /*!
     * Solve the linear system Ax=b by means of a Cholesky factorization.
     */
    void solve();
    //@}
}; // class EigenCholesky

//! Class EigenLU
/*!
 * @class EigenLU
 * This class implements a linear solver for the system Ax=b.
 * It uses the LU factorization on a square matrix.
 */
class EigenLU : public DirectSolver
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Empty constructor
    EigenLU(const UInt nbDofs = 0): 
		DirectSolver(nbDofs) {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
    EigenLU(const SpMat & A, const Vector & b):
        DirectSolver(A, b) {}
        
    //! Destructor
    ~EigenLU() = default;
    //@}

	//! @name Methods
    //@{
    //! Solve the linear system
    /*!
     * Solve the linear system Ax=b by means of a LU factorization.
     */
    void solve();
    //@}
}; // class EigenLU

//! Class EigenUmfPack
/*!
 * @class EigenUmfPack
 * This class implements a linear solver for the system Ax=b.
 * It uses the LU factorization on a square matrix (from UmfPack).
 */
class EigenUmfPack : public DirectSolver
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Empty constructor
    EigenUmfPack(const UInt nbDofs = 0): 
		DirectSolver(nbDofs) {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
    EigenUmfPack(const SpMat & A, const Vector & b):
        DirectSolver(A, b) {}
        
    //! Destructor
    ~EigenUmfPack() = default;
    //@}

    //! @name Methods
    //@{
    //! Solve the linear system
    /*!
     * Solve the linear system Ax=b by means of a UmfPack factorization.
     */
    void solve();
    //@}
};


//! Class IterativeSolver
/*!
 * @class IterativeSolver
 * This is an abstract class that implements a linear solver for the system Ax=b solved by means of an iterative method.
 * It stores the system matrix in a block form through a Saddle Point matrix, because only the blocks are needed to build
 * up the preconditioner and to solve the system.
 * Each derived class employs a different iterative solver.
 */
class IterativeSolver : public Solver
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Empty constructor
    IterativeSolver():
    Solver(), M_maxIter( S_referenceMaxIter ), M_iter(0), M_tol( S_referenceTol ), M_res(0), CIndex(0) {}
    
	//! Constructor
    /*!
     * @param Mdim M block dimension
     * @param Brow B block row
     * @param Bcol B block col
     */
	IterativeSolver(const UInt Mdim,const UInt Brow):
		M_A(Mdim,Brow), Solver(Mdim+Brow), 
		M_maxIter( S_referenceMaxIter ), M_iter(0), M_tol( S_referenceTol ), M_res(0), CIndex(0) {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
    IterativeSolver(const SaddlePointMat & A, const Vector & b ):
        M_A(A) ,Solver(b),
        M_maxIter( S_referenceMaxIter ), M_iter(0), M_tol( S_referenceTol ), M_res(0), CIndex(0) {}

    //! Destructor
    virtual ~IterativeSolver() = default;
    //@}

	//! @name Set Methods
    //@{
    //! Set the system matrix
    /*!
     * @param A the system matrix
     */
    void setA(const SaddlePointMat & A) { M_A = A; }
    
    //! Set the system dofs
    /*!
     * @param Mdim the M block dofs
     * @param Brow the B block rows
     * @param Bcol the B block cols
     */
    void setDofs(const UInt Mdim, const UInt Brow)
    {
		M_A.resize(Mdim, Brow);
		Solver::setDofs(Mdim+Brow);
    }
    
    //! Set the preconditioner
    /*!
     * @param precon the preconditioner type
     */
    void set_precon(const std::string prec)
    {
		preconPtr = preconHandler::Instance().getProduct(prec);
		preconPtr->set(M_A);  
    }
    
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

    //! @name Get Methods
    //@{
    //! Get the system matrix
    /*!
     * @return the system matrix
     */
    const SaddlePointMat & getA() const { return M_A; }

    //! Get the system matrix
    /*!
     * @return the system matrix
     */
    SaddlePointMat & getA() { return M_A; }
    
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
    
    //! Get the preconditioner
    /*!
     * @return a constant reference to the preconditioner
     */
    const preconditioner & getprec() const { return *preconPtr; }

    //! Get the preconditioner
    /*!
     * @return a reference to the preconditioner
     */
    preconditioner & getprec() { return *preconPtr; }

    //! Get the preconditioner pointer
    /*!
     * @return a pointer to the preconditioner
     */
    preconditioner * getprecPtr() { return preconPtr.get(); }
    
    //! Get the number of non zero
    /*!
     * @return the number of non zero
     */
    UInt getNonZero() { return M_A.nonZeros(); }
    //@}

    //! @name Methods
    //@{
    //! Print out the computation details
    /*!
     * Print out the cimputation details
     */
    virtual void print() const
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
    //@}

protected:
	//! The system matrix in Saddle Point form
	SaddlePointMat       M_A;
    //! Maximum number of iterations
    UInt                 M_maxIter;
    //! Number of iterations
    UInt                 M_iter;
    //! Tolerance
    Real                 M_tol;
    //! Residual error
    Real                 M_res;
    //! Computation index
    UInt                 CIndex;
    //! The preconditioner type
    preconPtr_Type       preconPtr;

    static constexpr Real S_referenceTol = 1e-6;
    static constexpr UInt S_referenceMaxIter = 100;
}; // class IterativeSolver

//! Class imlBiCGSTAB
/*!
 * @class imlBiCGSTAB
 * This class implements a linear solver for the system Ax=b.
 * It uses the stabilized bi-conjugate gradient method on a square matrix. There's also the possibility to use a restart 
 * for the 1st type of breakdown due to the fact that the actual residual becomes orthogonal to the initial one.
 * This strategy may lead to the convergence of the method, but it may also fail leading to a non convergent method
 * that does not stop itself due to 1st type breakdow.
 */
class imlBiCGSTAB : public IterativeSolver
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Empty constructor
    imlBiCGSTAB(): IterativeSolver(), restart(Default_restart) {}

	//! Constructor
    /*!
     * @param Mdim M block dimension
     * @param Brow B block row
     * @param Bcol B block col
     */
	imlBiCGSTAB(const UInt Mdim,const UInt Brow):
		IterativeSolver(Mdim,Brow), restart(Default_restart) {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
    imlBiCGSTAB( const SaddlePointMat & A, const Vector & b ):
        IterativeSolver(A, b), restart(Default_restart) {}
        
    //! Destructor
    ~imlBiCGSTAB() = default;
    //@}
    
    //! @name Set methods
    //@{
    //! Set the restart
    /*!
     * @param The restart
     */
    void setRestart(const bool rest){ restart = rest; };
    //@}
    
    //! @name Methods
    //@{
    //! Print out the computation details
    /*!
     * Print out the cimputation details
     */
    void print() const
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
    void solve();
    //@}
    
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
class imlFGMRES : public IterativeSolver
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Empty constructor
    imlFGMRES(): IterativeSolver(), m(Default_m) {}

	//! Constructor
    /*!
     * @param Mdim M block dimension
     * @param Brow B block row
     * @param Bcol B block col
     */
	imlFGMRES(const UInt Mdim,const UInt Brow):
		IterativeSolver(Mdim,Brow), m(Default_m) {}

    //! Constructor
    /*!
     * @param A Eigen sparse matrix
     * @param b RHS, it is Eigen vector
     */
    imlFGMRES( const SaddlePointMat & A, const Vector & b ):
        IterativeSolver(A, b), m(Default_m) {}
  
    //! Destructor
    ~imlFGMRES() = default;
    //@}
        
    //! @name Set Methods
    //@{  
    //! Set the restart
    /*!
     * @param The restart
     */
    void set_m(const UInt rest_val){ m = rest_val; };
	//@}
        
    //! @name Methods
    //@{
    //! Solve the linear system
    /*!
     * Solve the linear system Ax=b by means of the GMRES method.
     */
    void solve();
    //@}

private:
	//! The restart level (how many iterations needed to perform a restart)
	UInt m;
	//! The default restart value
	static constexpr UInt Default_m = 300;
	
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
    //! @name Constructor & Destructor
    //@{
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

    //! Destructor
    ~SamgSolver() = default;
	//@}

    //! @name Methods
    //@{
    //! Solve the linear system
   void solve();
   //@}

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
    //! @name Constructor & Destructor
    //@{
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
    //@}

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
    //! @name Constructor & Destructor
    //@{
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
    //@}

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
