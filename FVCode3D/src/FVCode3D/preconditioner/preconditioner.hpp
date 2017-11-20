/*!
 * @file preconditioner.hpp
 * @brief These classes implement simple classes that perform the "inversion" of the preconditioner.
 */

#ifndef __PRECONDITIONER_HPP__
#define __PRECONDITIONER_HPP__

#include <FVCode3D/core/BasicType.hpp>
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


//! Class for assembling a lumped inner product.
/*!
 * @class lump_InnerProduct
 * This class implements the lumping of the inner product matrix.
*/
class lump_InnerProduct
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Construct a diagonal preconditioner.
    /*!
     * @param Mat The matrix of the system that we want to precondition
     */
	lump_InnerProduct( const SpMat & ip_Mat):
		M(ip_Mat), M_lump(new DiagMat(ip_Mat.rows())) { (*M_lump).setZero(); }
	//! No Copy-Constructor
    lump_InnerProduct(const lump_InnerProduct &) = delete;
	//! No Empty-Constructor
    lump_InnerProduct() = delete;
	//! Destructor
    ~lump_InnerProduct() = default;
	//@}
	
	//! @name Get Methods
    //@{
	//! Get M lumped (read only)
    /*!
     * @return A const reference to the lumped inner product
     */
    const DiagMat & get_readOnly() const
        {return *M_lump;}
    
    //! Get M lumped 
    /*!
     * @return A reference to the lumped inner product
     */
    DiagMat & get() 
        {return *M_lump;}
    //@}
    
    //! @name Assemble Methods
    //@{
	//! Assemble the M lump
    /*!
     * Assemble the lumped inner product
     */
    void assemble();
    //@}    
    
private:
	//! A const reference to the inner product matrix
	const SpMat                & M;
	//! The diagonal matrix represented the lumped inner product
	std::unique_ptr<DiagMat>     M_lump;
};


//! Class for assembling a lumped inner product.
/*!
 * @class Inexact_SchurComplement
 * This class implements the inexact Schur Complement. It's a template class depending
 * on which type we choose for the inverse of the inner product matrix approximation.
*/
template<class MinvType>
class Inexact_SchurComplement
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Construct a diagonal preconditioner.
    /*!
     * @param Mat The matrix of the system that we want to precondition
     */
	Inexact_SchurComplement(const MinvType & Minverse, const SpMat & Bmat, const SpMat & Tmat):
		M_inv(Minverse), B(Bmat), T(Tmat), ISchurCompl(new SpMat(T.rows(),T.cols())) {}
	//! No Copy-Constructor
    Inexact_SchurComplement(const Inexact_SchurComplement &) = delete;
	//! No Empty-Constructor
    Inexact_SchurComplement() = delete;
	//! Destructor
    ~Inexact_SchurComplement() = default;
	//@}
	
	//! @name Get Methods
    //@{
	//! Get M lumped (read only)
    /*!
     * @return A const reference to the inexact Schur Complement
     */
    const SpMat & get_readOnly() const
        {return *ISchurCompl;}
    
    //! Get M lumped 
    /*!
     * @return A reference to the inexact Schur Complement
     */
    SpMat & get() 
        {return *ISchurCompl;}
    //@}
    
    //! @name Assemble Methods
    //@{
	//! Assemble the inexact Schur Complement
    /*!
     * Assemble the inexact Schur Complement
     */
    void assemble()
		{ *ISchurCompl = T - B*M_inv*B.transpose(); }
    //@}    
    
private:
	//! A const reference to the lumped inner product matrix
	const MinvType             & M_inv;
	//! A const reference to the divergence matrix
	const SpMat                & B;
	//! A const reference to the trasmissibility matrix
	const SpMat                & T;
	//! The inexact Schur Complement matrix
	std::unique_ptr<SpMat>       ISchurCompl;
};

}

#endif // __LOCAL_OPERATOR_HPP__
