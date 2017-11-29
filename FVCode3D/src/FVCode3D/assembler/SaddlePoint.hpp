/*!
 *  @file SaddlePoint.hpp
 *  @brief Class for building a saddle point problem.
 */

#ifndef __SADDLEPOINT_HPP__
#define __SADDLEPOINT_HPP__

#include <cmath>
#include <utility>
#include <FVCode3D/core/BasicType.hpp>
#include <FVCode3D/mesh/RigidMesh.hpp>
#include <FVCode3D/boundaryCondition/BC.hpp>

namespace FVCode3D
{

//! Base class for assembling a saddle point matrix
/*!
 * @class SaddlePoint_StiffMatrix
 * This is the class implementing the stiffness matrix as a block saddle point matrix,
 * useful in the case of MFD if we want to solve the system with an iterative
 * method preconditioning it in a not trivial way.
 */
class SaddlePoint_StiffMatrix
{
	
public:
	//! Typedef for std::pair<UInt,UInt>
	/*!
	* @typedef couple
	* This type definition permits to treat std::pair<UInt,UInt> as a couple.
	*/
	typedef std::pair<UInt,UInt> couple;
	
    //! @name Constructor & Destructor
    //@{
    //! Construct a SaddlePoint_Matrix, given a Rigid_Mesh.
    /*!
        @param rigid_mesh A Rigid_Mesh used to build the matrix
        @param size The size of the stiffness matrix
    */
    SaddlePoint_StiffMatrix(const Rigid_Mesh & pmesh, const BoundaryConditions & bc,
		const UInt mdim, const UInt brow, const UInt bcol, const UInt tdim):
        M_mesh(pmesh), M_bc(bc), Mdim(mdim), Bdim(brow,bcol), Tdim(tdim), 
        M(new SpMat(mdim,mdim)), B(new SpMat(brow,bcol)), T(new SpMat(tdim,tdim)),
        M_b(new Vector(Vector::Zero(mdim+brow))) {}
    //! No Copy-Constructor
    SaddlePoint_StiffMatrix(const SaddlePoint_StiffMatrix &) = delete;
    //! No Empty-Constructor
    SaddlePoint_StiffMatrix() = delete;
    //! Default Destructor
    virtual ~SaddlePoint_StiffMatrix() = default;
    //@}

    //! @name Get Methods
    //@{
    //! Get M block (read only)
    /*!
     * @return A const reference to the M block
     */
    const SpMat & getM_readOnly() const
        {return *M;}
    
    //! Get M block 
    /*!
     * @return A reference to the M block
     */
    SpMat & getM() 
        {return *M;}
        
    //! Get B block (read only)
    /*!
     * @return A const reference to the B block
     */
    const SpMat & getB_readOnly() const
        {return *B;}
    
    //! Get B block 
    /*!
     * @return A reference to the B block
     */
    SpMat & getB() 
        {return *B;}
        
    //! Get T block (read only)
    /*!
     * @return A const reference to the T block
     */
    const SpMat & getT_readOnly() const
        {return *T;}
    
    //! Get T block 
    /*!
     * @return A reference to the T block
     */
    SpMat & getT() 
        {return *T;}
        
    //! Get BC vector (read only)
    /*!
     * @return A reference to a constant vector that represents the part of the right hand side due to the boundary conditions.
     */
    const Vector & getBCVector_readOnly() const
        {return *M_b;}
        
    //! Get BC vector 
    /*!
     * @return A reference to a vector that represents the part of the right hand side due to the boundary conditions.
     */
    Vector & getBCVector() 
        {return *M_b;}

    //! Get size 
    /*!
     * @return The size of the whole matrix
     */
    UInt getSize() const
        {return Mdim+Bdim.first;}
        
    //! Get M size 
    /*!
     * @return The size of the M matrix
     */
    UInt getMSize() const
        {return Mdim;}
        
    //! Get size 
    /*!
     * @return The rows/column of B matrix
     */
    couple getBSize() const
        {return Bdim;}
    
    //! Get size 
    /*!
     * @return The size of the T matrix
     */
    UInt getTSize() const
        {return Tdim;}
    //@}

    //! @name Methods
    //@{
    //! Show system dimension 
    /*!
     * Show system dimension
     */
    void showMe() const
        {
			std::cout<<std::endl;
			std::cout<<"The system dimension is : "<<Mdim+Bdim.first<<std::endl<<std::endl;
		};
		
	//! Compress the block matrices
    /*!
     * Compress the block matrices M, B and T
     */
    void Compress() const
        {
			(*M).makeCompressed();
			(*B).makeCompressed();
			(*T).makeCompressed();
		};
        
    //! Assemble method
    /*!
     * Assemble the block matrices M, B and T
     */
    void assemble();
    //@}

private:
	//! Const reference to the rigid mesh
	const Rigid_Mesh              & M_mesh;
	//! Constant reference to the boundary conditions
	const BoundaryConditions      &	M_bc;  
	//! M block dimension
	UInt                            Mdim;
	//! B block row/column number
	couple                          Bdim;
	//! T block dimension    
	UInt                            Tdim;
    //! A unique ptr to the inner product matrix
    std::unique_ptr<SpMat>          M;
    //! A unique ptr to the Btilde block 
    std::unique_ptr<SpMat>          B;
    //! A unique pointer to Ttilde block
    std::unique_ptr<SpMat>          T;
	//! A reference to the rhs of the system
	std::unique_ptr<Vector>         M_b;
};

} //FVCode3D

#endif //__SADDLEPOINT_HPP__

