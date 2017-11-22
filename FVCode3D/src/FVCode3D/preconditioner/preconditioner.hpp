/*!
 * @file preconditioner.hpp
 * @brief These classes implement simple classes that perform the "inversion" of the preconditioner.
 */

#ifndef __PRECONDITIONER_HPP__
#define __PRECONDITIONER_HPP__

#include <FVCode3D/core/BasicType.hpp>
#include <FVCode3D/assembler/SaddlePoint.hpp>
#include <Eigen/LU>
#include <unsupported/Eigen/SparseExtra>

namespace FVCode3D
{
	
//! Base class for assembling a preconditioner.
/*!
 * @class preconditioner
 * This is the base class for a preconditioner. As derived classes we can have any preconditioner:
 * identity, diagonal, inexact LU factorization and so on. The goal of these classes is not the
 * assemblation of the preconditioner matrix, but solving the linear system Pz = r in the iterative
 * scheme. They act as a wrapper around that allow the practical usage of preconditioner.
*/
class preconditioner
{
public:
    //! @name Solve Methods
    //@{
    //! Solve the linear system Pz=r.
    /*!
     * @param r The rhs vector on what we apply the P^-1
     * Note that the return value is passed through move semantic.
     */
    virtual Vector solve(const Vector & r) const = 0;
    //@}                  	
};


//! Class for assembling an identity preconditioner
/*!
 * @class identity_preconditioner
 * This class constructs an identity preconditioner. Create and pass an object of this type
 * in the case of running an iterative solver without preconditioning.
 */
class identity_preconditioner: public preconditioner
{
public:
    //! @name Solve Methods
    //@{
    //! Solve the linear system Pz=r
    /*!
     * @param r The rhs vector on what we apply the P^-1
     */
    Vector solve(const Vector & r) const
		{ return r; };
    //@}                 
};


//! Class for assembling a diagonal preconditioner
/*!
 * @class diagonal_preconditioner
 * This class constructs a diagonal preconditioner.
 */
class diagonal_preconditioner: public preconditioner
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Construct a diagonal preconditioner.
    /*!
     * @param Mat The matrix of the system that we want to precondition
     */
	diagonal_preconditioner( const SpMat & Mat ):
		A(Mat){}
	//! No Copy-Constructor
    diagonal_preconditioner(const diagonal_preconditioner &) = delete;
	//! No Empty-Constructor
    diagonal_preconditioner() = delete;
	//! Destructor
    ~diagonal_preconditioner() = default;
	//@}

    //! @name Solve Methods
    //@{
    //! Solve the linear system Pz=r
    /*!
     * @param r The rhs vector on what we apply the P^-1
     */
    Vector solve(const Vector & r) const;
    //@}
    
private:            
	//! The matrix of the system that we want to precondition
    const SpMat &   A;                  
};


//! Class for assembling a lumped inner product builder.
/*!
 * @class lumpIP_builder
 * This class implements the lumping of the inner product matrix.
*/
class lumpIP_builder
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Construct a lumped inner product builder.
    /*!
     * @param ip_Mat The innner product matrix
     */
	lumpIP_builder( const SpMat & ip_Mat):
		M(ip_Mat) {}
	//! No Copy-Constructor
    lumpIP_builder(const lumpIP_builder &) = delete;
	//! No Empty-Constructor
    lumpIP_builder() = delete;
	//! Destructor
    ~lumpIP_builder() = default;
	//@}
    
    //! @name Assemble Methods
    //@{
	//! Assemble the lumped inner product
    /*!
     * @param M_lump A reference to the lumped inner product to be built
     */
    void build(DiagMat & M_lump) const;
    //@}    
    
private:
	//! A const reference to the inner product matrix
	const SpMat       & M;
};


//! Class for assembling an inexact Schur Complement builder.
/*!
 * @class ISchurComplement_builder
 * This class implements the inexact Schur Complement with a lumping to approximate the 
 * inverse of the inner product matrix.
*/
class ISchurComplement_builder
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Construct an inexact Schur Complement.
    /*!
     * @param M_lump The lumped inner product matrix
     * @param Bmat The B block of the saddle point system
     * @param Tmat The T block of the saddle point system
     */
	ISchurComplement_builder(const DiagMat & Ml, const SpMat & Bmat, const SpMat & Tmat):
		Mlump(Ml), B(Bmat), T(Tmat) {}
	//! No Copy-Constructor
    ISchurComplement_builder(const ISchurComplement_builder &) = delete;
	//! No Empty-Constructor
    ISchurComplement_builder() = delete;
	//! Destructor
    ~ISchurComplement_builder() = default;
	//@}
    
    //! @name Assemble Methods
    //@{
	//! Assemble the inexact Schur Complement
    /*!
     * @param A reference to the inexact Schur Complement to be built
     */
    void build(SpMat & ISchurCompl) const
		{
			ISchurCompl =  - B * Mlump.inverse() * B.transpose();
			ISchurCompl += T; 
		}
    //@}    
    
private:
	//! A const reference to the lumped inner product matrix
	const DiagMat              & Mlump;
	//! A const reference to the B block 
	const SpMat                & B;
	//! A const reference to the T block
	const SpMat                & T;
};


//! Class for assembling a saddle point matrix.
/*!
 * @class SPMatrix
 * This class implements a generic saddle point matrix. It consists of 3 block matrices
 * M, B and T. It is build up with through an object of type SaddlePoint_StiffMatrix
 * that is the assembler of the numerical method and it overloads the *(Vector) operator
 * to use the blocks matrix to compute the matrix-vector product. In this way we avoid
 * to storage the whole system matrix and we use the same 3 blocks (without any copy of them)
 * to make everything we need: matrix-vector product, building and inverting the preconditioner.
 * Clearly this class is interesting only with an iterative system solving the system.
 * Constructors, assignement, destructor and so on are all defaulted.
*/
class SPMatrix
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Construct a saddle point matrix
    /*!
     * @param SP_Stiff A reference to the saddle point stiffness matrix
     */
	SPMatrix(SaddlePoint_StiffMatrix & SP_Stiff):
		M(SP_Stiff.getM()), B(SP_Stiff.getB()), T(SP_Stiff.getT()) {}
		//! No Copy-Constructor
    SPMatrix(const SPMatrix &) = default;
	//! No Empty-Constructor
    SPMatrix() = default;
	//! Destructor
    ~SPMatrix() = default;
	//@}
		
    //! @name Set Methods
    //@{
	//! Set the saddle point matrix
    /*!
     * @param SP_Stiff A reference to the saddle point stiffness matrix
     */
    void Set(SaddlePoint_StiffMatrix & SP_Stiff)
	{
		M = SP_Stiff.getM();
		B = SP_Stiff.getB();
		T = SP_Stiff.getT();
	}
    //@}
    
    //! @name Get Methods
    //@{
    //! Get M block (read only)
    /*!
     * @return A const reference to the M block
     */
    const SpMat & getM() const
        {return M;}
    
    //! Get B block 
    /*!
     * @return A reference to the B block
     */
    const SpMat & getB() const
        {return B;}
        
    //! Get T block (read only)
    /*!
     * @return A const reference to the T block
     */
    const SpMat & getT() const
        {return T;}
    //@}

    //! @name Operators
    //@{
	//! Overload of matrix-vector product operator using the blocks
    /*!
     * @param x A const reference to the vector
     * @return The vector resulting from the matrix-vector product 
     */
    Vector operator * (const Vector & x) const
    {
		Vector result(M.rows()+B.rows());
		result.segment(0, M.rows()) = M*x.segment(0,M.rows()) + B.transpose()*x.segment(M.rows(),B.rows());
		result.segment(M.rows(),B.rows()) = B*x.segment(0,M.rows()) + T*x.segment(M.rows(),B.rows());
		return result;
	}    
    //@}

private:
	//! The M block matrix
	SpMat      M;
	//! The B block matrix 
	SpMat      B;
	//! The T block matrix
	SpMat      T;
};


//! Class for assembling a saddle point matrix.
/*!
 * @class SPMatrix
 * This class builds up a block triangolar preconditioner where the M and SC
 * blocks are approximated using the lumping on the inner product matrix.
 * The SC linear system is solved through Cholesky.
*/
class BlockTriangular_preconditioner: public preconditioner
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Construct a diagonal preconditioner.
    /*!
     * @param SP The saddle point matrix
     */
	BlockTriangular_preconditioner( const SPMatrix & SP ): B(SP.getB()){}
	//! No Copy-Constructor
    BlockTriangular_preconditioner(const BlockTriangular_preconditioner &) = delete;
	//! No Empty-Constructor
    BlockTriangular_preconditioner() = delete;
	//! Destructor
    ~BlockTriangular_preconditioner() = default;
	//@}

    //! @name Assemble Methods
    //@{
	//! Assemble the approximations of M and SC
    /*!
     * @param SP_Stiff A reference to the saddle point stiffness matrix
     */
    void assemble(const SPMatrix & SP)
	{
		lumpIP_builder lumpBuild(SP.getM());
		lumpBuild.build(IM);
		ISchurComplement_builder ISCBuilder(IM, B, SP.getT());
		ISCBuilder.build(ISC);
	}
    //@}

    //! @name Solve Methods
    //@{
    //! Solve the linear system Pz=r
    /*!
     * @param r The rhs vector on what we apply the P^-1
     */
    Vector solve(const Vector & r) const;
    //@}
    
private:
	//! The B block matrix
    const SpMat         & B;
	//! The M approximation
    DiagMat               IM;   
	//! The inexact Schur Complement matrix
    SpMat                 ISC;                           
};


}

#endif // __LOCAL_OPERATOR_HPP__
