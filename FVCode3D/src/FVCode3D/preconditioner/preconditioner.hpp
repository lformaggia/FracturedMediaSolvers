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
	
}

#endif // __LOCAL_OPERATOR_HPP__
