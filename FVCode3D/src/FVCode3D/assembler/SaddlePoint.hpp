/*!
 *  @file SaddlePoint.hpp
 *  @brief Class for building a saddle point problem.
 */

#ifndef __SADDLEPOINT_HPP__
#define __SADDLEPOINT_HPP__

#include <cmath>
#include <utility>
#include <unsupported/Eigen/SparseExtra>
#include <FVCode3D/core/BasicType.hpp>
#include <FVCode3D/mesh/RigidMesh.hpp>
#include <FVCode3D/preconditioner/preconditioner.hpp>
#include <FVCode3D/boundaryCondition/BC.hpp>

namespace FVCode3D
{

//! Base class for assembling a saddle point stiffness matrix
/*!
 * @class SaddlePoint_StiffMatHandler
 * This is the class that catually build the stiffness matrix as a block saddle point matrix,
 * useful in the case of MFD if we want to solve the system with an iterative
 * method preconditioning it in a not trivial way.
 */
class SaddlePoint_StiffMatHandler
{
	
public:
	
    //! @name Constructor & Destructor
    //@{
    //! Construct a SaddlePoint_StiffMatHandler.
    /*!
        @param rigid_mesh A Rigid_Mesh used to build the matrix
        @param bc A boundary condition
        @param Msp The saddle point matrix
        @param b The rhs
    */
    SaddlePoint_StiffMatHandler(const Rigid_Mesh & pmesh, SaddlePointMat & Msp, Vector & b, const BoundaryConditions & bc):
        M_mesh(pmesh), M_SP(Msp), M_b(b), M_bc(bc){}
    //! No Copy-Constructor
    SaddlePoint_StiffMatHandler(const SaddlePoint_StiffMatHandler &) = delete;
    //! No Empty-Constructor
    SaddlePoint_StiffMatHandler() = delete;
    //! Default Destructor
    virtual ~SaddlePoint_StiffMatHandler() = default;
    //@}

    //! @name Get Methods
    //@{
    //! Get M block (read only)
    /*!
     * @return A const reference to the M block
     */
    const SpMat & getM() const
        {return M_SP.getM();}
    
    //! Get M block 
    /*!
     * @return A reference to the M block
     */
    SpMat & getM() 
        {return M_SP.getM();}
        
    //! Get B block (read only)
    /*!
     * @return A const reference to the B block
     */
    const SpMat & getB() const
        {return M_SP.getB();}
    
    //! Get B block 
    /*!
     * @return A reference to the B block
     */
    SpMat & getB() 
        {return M_SP.getB();}
        
    //! Get T block (read only)
    /*!
     * @return A const reference to the T block
     */
    const SpMat & getT() const
        {return M_SP.getT();}
    
    //! Get T block 
    /*!
     * @return A reference to the T block
     */
    SpMat & getT() 
        {return M_SP.getT();}
        
    //! Get BC vector (read only)
    /*!
     * @return A reference to a constant vector that represents the part of the right hand side due to the boundary conditions.
     */
    const Vector & getBC() const
        {return M_b;}
        
    //! Get BC vector 
    /*!
     * @return A reference to a vector that represents the part of the right hand side due to the boundary conditions.
     */
    Vector & getBC() 
        {return M_b;}

    //! Get size 
    /*!
     * @return The size of the whole matrix
     */
    UInt getSize() const
        {return M_SP.getM().rows()+M_SP.getB().rows();}
        
    //! Get M size 
    /*!
     * @return The size of the M matrix
     */
    UInt getMsize() const
        {return M_SP.getM().rows();}
        
    //! Get B rows
    /*!
     * @return The rows of B matrix
     */
    UInt getBrows() const
        {return M_SP.getB().rows();}
    
    //! Get B cols
    /*!
     * @return The column of B matrix
     */
    UInt getBcols() const
        {return M_SP.getM().rows();}
    
    //! Get T size 
    /*!
     * @return The size of the T matrix
     */
    UInt getTsize() const
        {return M_SP.getB().rows();}
    //@}

    //! @name Methods
    //@{
    //! Export method 
    /*!
     * export the M block
     */
    void ExportM()
    {
		Eigen::saveMarket( M_SP.getM(), "Mblock.m" );
	}
	
	//! Show method 
    /*!
     * show the system dimension
     */
    void show()
    {
		std::cout<<std::endl;
		std::cout<<"The system dimension is : "<<getSize()<<std::endl<<std::endl;
	}
	    
	//! Set dofs 
    /*!
     * @param size The dofs to be set
     */
    void setDofs(const UInt Mdim, const UInt Brow)
    {
		M_SP.resize(Mdim,Brow);
		M_b.resize(Mdim+Brow);
	}
        
    //! Assemble method
    /*!
     * Assemble the block matrices M, B and T
     */
    void assemble();
    //@}

private:
	//! Const reference to the rigid mesh
	const Rigid_Mesh              & M_mesh; 
    //! A reference to the saddle point matrix
    SaddlePointMat                & M_SP;
	//! A reference to the rhs of the system
	Vector                        & M_b;
	//! Constant reference to the boundary conditions
	const BoundaryConditions      &	M_bc; 
};

} //FVCode3D

#endif //__SADDLEPOINT_HPP__

