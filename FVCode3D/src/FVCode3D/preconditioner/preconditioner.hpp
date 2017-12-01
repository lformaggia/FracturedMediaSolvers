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
 * scheme. They act as a wrapper that allows the practical usage of preconditioner.
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
 * This class implements the inexact Schur Complement with the inverse of the diagonal block M
 * to approximate the inverse of the inner product matrix.
*/
class ISchurComplement_builder
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Construct an inexact Schur Complement.
    /*!
     * @param Mdi The diagonal inverse matrix
     * @param Bmat The B block of the saddle point system
     * @param Tmat The T block of the saddle point system
     */
	ISchurComplement_builder(const DiagMat & Mdi, const SpMat & Bmat, const SpMat & Tmat):
		Md_inv(Mdi), B(Bmat), T(Tmat) {}
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
			ISchurCompl =  - B * Md_inv * B.transpose();
			ISchurCompl += T; 
		}
    //@}    
    
private:
	//! A const reference to the inner product matrix
	const DiagMat              & Md_inv;
	//! A const reference to the B block 
	const SpMat                & B;
	//! A const reference to the T block
	const SpMat                & T;
};


//! Class for assembling a saddle point matrix.
/*!
 * @class SaddlePointMat
 * This class implements a generic saddle point matrix. It consists of 3 block matrices
 * M, B and T. It is build up through an object of type SaddlePoint_StiffMatrix
 * that is the assembler of the numerical method and it overloads the *(Vector) operator
 * to use the blocks matrix to compute the matrix-vector product. In this way we avoid
 * to storage the whole system matrix and we use the same 3 blocks (without any copy of them)
 * to make everything we need: matrix-vector product, building and inverting the preconditioner.
 * Clearly this class is interesting only with an iterative system solving the system.
 * Constructors, assignement, destructor and so on are all defaulted.
*/
class SaddlePointMat
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Construct a saddle point matrix
    /*!
     * @param SP_Stiff A reference to the saddle point stiffness matrix
     */
	SaddlePointMat(SaddlePoint_StiffMatrix & SP_Stiff):
		M(SP_Stiff.getM()), B(SP_Stiff.getB()), T(SP_Stiff.getT()) {}
	//! Copy-Constructor
    SaddlePointMat(const SaddlePointMat &) = default;
	//! Empty-Constructor
    SaddlePointMat() = default;
	//! Destructor
    ~SaddlePointMat() = default;
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
		result.segment(0,M.rows()) = M*x.segment(0,M.cols()) + B.transpose()*x.segment(M.cols(),B.rows());
		result.segment(M.rows(),B.rows()) = B*x.segment(0,M.cols()) + T*x.segment(M.cols(),B.rows());
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
	diagonal_preconditioner( const SaddlePointMat & Mat ):
		M(Mat.getM()), T(Mat.getT()) {}
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
	//! A const reference to the inner product matrix
    const SpMat &       M;    
	//! A const reference to the T block matrix
    const SpMat &       T;               
};

//! Base class for assembling a Block Triangular preconditioner.
/*!
 * @class BlockTriangular_preconditioner
 * Class that builds up a block triangolar preconditioner. The M block is
 * approximated through the diagonal part of M and the SC through the inverse
 * of the diagonal part of M. The SC system is solved through the CG method.
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
	BlockTriangular_preconditioner( const SaddlePointMat & SP ): 
		B(SP.getB()), Md_inv(SP.getM().rows()), ISC(SP.getB().rows(),SP.getB().rows()), MaxIt(MaxIt_Default), tol(tol_Default) {}
	//! No Copy-Constructor
    BlockTriangular_preconditioner(const BlockTriangular_preconditioner &) = delete;
	//! No Empty-Constructor
    BlockTriangular_preconditioner() = delete;
	//! Destructor
    ~BlockTriangular_preconditioner() = default;
	//@}
	
	//! @name Get Methods
    //@{
    //! Get ISC block 
    /*!
     * @return A reference to the ISC block
     */
	const SpMat & getISC() const
		{return ISC;}
    //@}
    
    //! @name Set Methods
    //@{
    //! Set the max it value for CG
    /*!
     * @param itmax Max it value for CG
     */
	void setMaxIt(const UInt itmax)
		{MaxIt = itmax;}
		
	//! Set the tolerance value for CG
    /*!
     * @param Tol tolerance value for CG
     */
	void set_tol(const UInt Tol)
		{tol = Tol;}
    //@}

    //! @name Assemble Methods
    //@{
	//! Assemble the the inverse diag of M and the SC
    /*!
     * @param SP_Stiff A reference to the saddle point stiffness matrix
     */
    void assemble(const SaddlePointMat & SP)
	{
		Md_inv = SP.getM().diagonal().asDiagonal().inverse();
		ISchurComplement_builder ISCBuilder(Md_inv, B, SP.getT());
		ISCBuilder.build(ISC);
	}
    //@}

    //! @name Solve Methods
    //@{
    //! Solve the linear system Pz=r
    /*!
     * @param r The rhs vector on what we apply the P^-1
     */
    virtual Vector solve(const Vector & r) const;
    //@}
    
protected:
	//! The B block matrix
    const SpMat              & B;
    //! The inverse of the diagonal of M
    DiagMat                    Md_inv;
	//! The inexact Schur Complement matrix
    SpMat                      ISC;
    //! The max it for CG
    UInt                       MaxIt;
    //! The tolerance for CG
    Real                       tol;
    //! The max it for CG (default value)
    static constexpr UInt      MaxIt_Default = 200;
    //! The tolerance for CG (default value)
    static constexpr Real      tol_Default = 1e-2;
                      
};

//! Class for assembling a ILU preconditioner.
/*!
 * @class ILU_preconditioner
 * This class builds up a ILU preconditioner where the M block is approximated via its diagonal part
 * and the SC is approximated through the inverse of the diagonal part of M.
 * The SC linear system is solved through the Eigen CG.
*/
class ILU_preconditioner: public preconditioner
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Construct a diagonal preconditioner.
    /*!
     * @param SP The saddle point matrix
     */
	ILU_preconditioner( const SaddlePointMat & SP ): 
		B(SP.getB()), Md_inv(SP.getM().rows()), ISC(SP.getB().rows(),SP.getB().rows()), MaxIt(MaxIt_Default), tol(tol_Default)
			{ Md_inv.setZero(); }
	//! No Copy-Constructor
    ILU_preconditioner(const ILU_preconditioner &) = delete;
	//! No Empty-Constructor
    ILU_preconditioner() = delete;
	//! Destructor
    ~ILU_preconditioner() = default;
	//@}
	
	//! @name Get Methods
    //@{
    //! Get B block 
    /*!
     * @return A reference to the B block
     */
	const SpMat & getISC() const
		{return ISC;}
    //@}
    
    //! @name Set Methods
    //@{
    //! Set the max it value for CG
    /*!
     * @param itmax Max it value for CG
     */
	void setMaxIt(const UInt itmax)
		{MaxIt = itmax;}
		
	//! Set the tolerance value for CG
    /*!
     * @param Tol tolerance value for CG
     */
	void set_tol(const UInt Tol)
		{tol = Tol;}
    //@}

    //! @name Assemble Methods
    //@{
	//! Assemble the approximations of M and SC
    /*!
     * @param SP_Stiff A reference to the saddle point stiffness matrix
     */
    void assemble(const SaddlePointMat & SP)
	{
		Md_inv = SP.getM().diagonal().asDiagonal().inverse();
		ISchurComplement_builder ISCBuilder(Md_inv, B, SP.getT());
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
    const SpMat              & B;
    //! The inverse of the diagonal of M
    DiagMat                    Md_inv; 
	//! The inexact Schur Complement matrix
    SpMat                      ISC;
    //! The max it for CG
    UInt                       MaxIt;
    //! The tolerance for CG
    Real                       tol;
    //! The max it for CG (default value)
    static constexpr UInt      MaxIt_Default = 200;
    //! The tolerance for CG (default value)
    static constexpr Real      tol_Default = 1e-2;
                      
};

}

#endif // __LOCAL_OPERATOR_HPP__
