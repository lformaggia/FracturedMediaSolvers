/*!
 *  @file MatrixHandler.hpp
 *  @brief Base Class for building a matrix by discretizing a PDE with a finite volume method.
 */

#ifndef __MATRIXHANDLER_HPP__
#define __MATRIXHANDLER_HPP__

#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
#include <FVCode3D/core/TypeDefinition.hpp>
#include <FVCode3D/mesh/RigidMesh.hpp>
#include <FVCode3D/core/Data.hpp>

namespace FVCode3D
{

//! Base class for assembling a system matrix
/*!
 * @class MatrixHandler
 * This class is a base class used as a model for the derived classes (such as the FV stiffness matrix or the MFD one).
 * It implements a square-Matrix (N x N). It needs a Rigid_Mesh and a size to be built.
 * Abstract class.
 */
class MatrixHandler
{
	
public:
    //! @name Constructor & Destructor
    //@{

    //! Construct a MatrixHandler, given a Rigid_Mesh.
    /*!
        @param rigid_mesh A Rigid_Mesh used to build the matrix
        @param size The size of the stiffness matrix
    */
    MatrixHandler(const Rigid_Mesh & rigid_mesh, SpMat & Mat):
        M_mesh(rigid_mesh), M_Matrix(Mat){}
			
    //! No Copy-Constructor
    MatrixHandler(const MatrixHandler &) = delete;
    //! No Empty-Constructor
    MatrixHandler() = delete;
    //! Default Destructor
    virtual ~MatrixHandler() = default;
    //@}

    //! @name Get Methods
    //@{
    //! Get Matrix (read only)
    /*!
     * @return A const reference to the matrix
     */
    const SpMat & getMatrix() const
        {return M_Matrix;}
    
    //! Get Matrix 
    /*!
     * @return A reference to the matrix
     */
    SpMat & getMatrix() 
        {return M_Matrix;}

    //! Get size (const)
    /*!
     * @return The size N (number of rows) of the matrix
     */
    virtual UInt getSize() const
        {return M_Matrix.size();}
    //@}

    //! @name Methods
    //@{
    //! Set dofs 
    /*!
     * @param size The dofs to be set
     */
    virtual void setDofs(const UInt size)
    {
		M_Matrix.resize(size,size);
	}
        
    //! Assemble method
    /*!
     * Assemble the matrix
     */
    virtual void assemble()=0;
    //@}

protected:
    //! A constant reference to the Rigid_Mesh
    const Rigid_Mesh & M_mesh;
    //! A reference to the system matrix
    SpMat & M_Matrix;
};


//! Base class for assembling a FV system matrix
/*!
 * @class MatrixHandlerFV
 * This class is a base class used as a model for the derived 
 * FV classes (such as the FV stiffness matrix or FV mass matrix).
 * The used scheme is based on the Eigen triplets.
 * Abstract class.
 */
class MatrixHandlerFV: public MatrixHandler
{
	public:
    //! @name Constructor & Destructor
    //@{
    //! Construct a MatrixHandler, given a Rigid_Mesh.
    /*!
        @param rigid_mesh A Rigid_Mesh used to build the matrix
        @param M_numet It is the numerical method used for the discretization
    */
    MatrixHandlerFV(const Rigid_Mesh & rigid_mesh, SpMat & Mat):
        MatrixHandler(rigid_mesh, Mat), M_offsetRow(0), M_offsetCol(0){}
    //! No Copy-Constructor
    MatrixHandlerFV(const MatrixHandlerFV &) = delete;
    //! No Empty-Constructor
    MatrixHandlerFV() = delete;
    //! Default Destructor
    virtual ~MatrixHandlerFV() = default;
    //@}

    //! @name Get Methods
    //@{
    //! Get Triplets (const)
    /*!
     * @return A reference to the triplets
     */
    const std::vector<Triplet> & getTriplets() const
        {return M_matrixElements;}

    //! Get size (const)
    /*!
     * @return The size N (number of rows) of the matrix
     */
    UInt getSize() const
        {return M_Matrix.size() + M_offsetRow;}
    //@}

    //! @name Assemble Methods
    //@{
    //! Assemble method
    /*!
     * Assemble the matrix
     */
    virtual void assemble()=0;
    //@}

    //! @name Methods
    //@{
    //! Close the matrix
    /*!
     * Fill the matrix with triplets and clear the vector of triplets
     */
    virtual void closeMatrix()
    {
        M_Matrix.setFromTriplets( std::begin( M_matrixElements ), std::end( M_matrixElements ) );
        M_matrixElements.clear();
    }

    //! Clear triplets
    /*!
     * Clear the vector of triplets
     */
    void clearTriplets()
    {
        M_matrixElements.clear();
    }

    //! Set offsets
    /*!
     * Set offsets and resize the matrix as number of dofs + offsets
     * @param row row offset
     * @param col column offset
     */
    virtual void setOffsets(const UInt row, const UInt col)
    {
        M_offsetRow = row;
        M_offsetCol = col;
        M_Matrix = SpMat(M_Matrix.size() + M_offsetRow, M_Matrix.size() + M_offsetCol);
    }
    //@}

protected:
    //! Row offset
    UInt M_offsetRow;
    //! Column offset
    UInt M_offsetCol;
    //! Vector of triplets
    std::vector<Triplet> M_matrixElements;
};

//! Base class for assembling a MFD system matrix
/*!
 * @class MatrixHandlerMFD
 * This class is a base class used as a model for the derived system 
 * MFD classes (such as the MFD stiffness matrix).
 * We build up the matrix system directly with insert, so we do not need the
 * triplets but we need to compress it after having built it.
 * Abstract class.
 */
class MatrixHandlerMFD: public MatrixHandler
{
	public:
    //! @name Constructor & Destructor
    //@{
    //! Construct a MatrixHandler, given a Rigid_Mesh.
    /*!
        @param rigid_mesh A Rigid_Mesh used to build the matrix
        @param M_numet It is the numerical method used for the discretization
    */
    MatrixHandlerMFD(const Rigid_Mesh & rigid_mesh, SpMat & Mat):
		MatrixHandler(rigid_mesh, Mat){}
    //! No Copy-Constructor
    MatrixHandlerMFD(const MatrixHandlerMFD &) = delete;
    //! No Empty-Constructor
    MatrixHandlerMFD() = delete;
    //! Default Destructor
    virtual ~MatrixHandlerMFD() = default;
    //@}

    //! @name Assemble Methods
    //@{
    //! Assemble method
    /*!
     * Assemble the matrix
     */
    virtual void assemble()=0;
    //@}
};

} // namespace FVCode3D

#endif // __MATRIXHANDLER_HPP__
