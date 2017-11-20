/*!
 * @file Stiffness.hpp
 * @brief This class build a Stiffness-matrix of the Darcy problem.
 */

#ifndef __DARCYSTIFFNESS_HPP__
#define __DARCYSTIFFNESS_HPP__

#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
#include <FVCode3D/assembler/MatrixHandler.hpp>
#include <FVCode3D/boundaryCondition/BC.hpp>
#include <FVCode3D/core/Data.hpp>

namespace FVCode3D
{

class Rigid_mesh;
class PropertiesMap;

//! Class for assembling a stiffness matrix
/*!
 * @class StiffMatrixFV
 * This class constructs the FV stiffness-matrix for the Darcy problem.
 * The adopted technique is a two point finite volume method.
 * The fractures are considered as cells and take part to discretization.
 */
class StiffMatrixFV: public MatrixHandlerFV
{

    //! Typedef for std::pair<UInt,UInt>
    /*!
     * @typedef Fracture_Juncture
     * This type definition permits to treat std::pair<UInt,UInt> as a Fracture_Juncture.
     */
    typedef std::pair<UInt,UInt> Fracture_Juncture;
    //! Typedef for Rigid_Mesh::Facet_ID
    /*!
     * @typedef Facet_ID
     * This type definition permits to treat Rigid_Mesh::Facet_ID as a Facet_ID.
     */
    typedef Rigid_Mesh::Facet_ID Facet_ID;
    //! Typedef for Rigid_Mesh::Edge_ID
    /*!
     * @typedef Edge_ID
     * This type definition permits to treat Rigid_Mesh::Edge_ID as a Edge_ID.
     */
    typedef Rigid_Mesh::Edge_ID Edge_ID;

public:
    //! @name Constructor & Destructor
    //@{

    //! Construct a stiffness-Matrix, given a Rigid_Mesh and the boundary conditions
    /*!
     * @param rigid_mesh A Rigid_Mesh used to build the matrix
     * @param BC Boundary conditions given in the container BoundaryConditions
     */
    StiffMatrixFV(const Rigid_Mesh & rigid_mesh, UInt size, const BoundaryConditions & BC):
        MatrixHandlerFV(rigid_mesh, size), M_b(new Vector(Vector::Zero(size))),
        M_bc(BC) {}
    //! No Copy-Constructor
    StiffMatrixFV(const StiffMatrixFV &) = delete;
    //! No Empty-Constructor
    StiffMatrixFV() = delete;
    //! Destructor
    ~StiffMatrixFV() = default;
    //@}

    //! @name Get Methods
    //@{
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
    //@}

    //! @name Methods
    //@{
    //! Assemble method
    /*!
     * Assemble the stiffness matrix
     */
    void assemble();

    //! Set offsets
    /*!
     * Set offsets and resize the matrix as number of dofs + offsets.
     * Further, it resizes the RHS and initializes it to zero.
     * @param row row offset
     * @param col column offset
     */
    virtual void setOffsets(const UInt row, const UInt col)
    {
        this->setOffsets(row, col);
        M_b.reset( new Vector( Vector::Constant( this->M_size + this->M_offsetRow, 0. ) ) );
    }

    //@}

public:

    //! @name Alpha Methods
    //@{

    //! Border center
    /*!
     * @param fj the Id of the juncture of two Fracture_Facet in 3D
     * @return The center of the juncture between two Fracure_Facet
     */
    Point3D getBorderCenter(Fracture_Juncture fj) const;

    //! It is called by the method assemble() and it computes the coefficient alpha
    /*!
     * @param cellId the Id of a Cell
     * @param facet A pointer to a Rigid_mesh::Facet_ID
     * @return The computed coefficient alpha
     */
    Real findAlpha (const UInt & cellId, const Facet_ID * facet) const;

    //! It is called by the method assemble() and it computes the coefficient alpha
    /*!
     * @param facetId the Id of a Facet
     * @param edge A pointer to a Rigid_mesh::Edge_ID
     * @return The computed coefficient alpha
     */
    Real findAlpha (const UInt & facetId, const Edge_ID * edge) const;

    //! It is called by the method assemble() and it computes the coefficient alpha in the case of Dirichlet BC
    /*!
     * @param cellId the Id of a Cell
     * @param facet A pointer to a Rigid_mesh::Facet_ID
     * @return The computed coefficient alpha
     */
    Real findDirichletAlpha (const UInt & cellId, const Facet_ID * facet) const;

    //! It is called by the method assemble() and it computes the coefficient alpha in the case of Dirichlet BC
    /*!
     * @param facetId the Id of a Facet
     * @param edge A pointer to a Rigid_mesh::Edge_ID
     * @return The computed coefficient alpha
     */
    Real findDirichletAlpha (const UInt & facetId, const Edge_ID * edge) const;

    //! It is called by the method assemble() and it computes the coefficient alpha in the case of a fracture in 3D
    /*!
     * @param fj is a Fracture_Juncture
     * @param n_Id The Id of the Fracture_Facet
     * @return The computed coefficient alpha
     */
    Real findFracturesAlpha (const std::pair<UInt,UInt> fj, const UInt n_Id) const;
    //@}

protected:

    //! @name Assemble Methods
    //@{

    //! Assemble the porous medium block
    /*!
     * Assemble the porous media block
     */
    void assemblePorousMatrix();

    //! Assemble the BCs for the porous medium
    /*!
     * Assemble the BCs for the porous medium
     */
    void assemblePorousMatrixBC();

    //! Assemble the fractures block
    /*!
     * Assemble the fractures block
     */
    void assembleFracture();

    //! Assemble the BCs for the fractures
    /*!
     * Assemble the BCs for the fractures
     */
    void assembleFractureBC();
    //@}

protected:
    //! Unique pointer to the vector that contains the effects of BCs on the RHS
    std::unique_ptr<Vector> M_b;
    //! The constant container of the BCs
    const BoundaryConditions & M_bc;
};


//! Class for assembling a MFD system matrix.
/*!
 * @class StiffMatrixMFD
 * This is the class to assemble the MFD stiffness matrix. It builds up the matrix system S
 * in an efficient way using the global_BulkBuilder and the FractureBuilder. The bulk builder
 * builds up the M, B, Dt matrices, then we impose the bulk bcs on M and on the rhs. The 
 * fracture builder builds up the coupling conditions matrices C and Ct, it modifies properly
 * the M matrix (coupling conditions again) and builds up the trasmissibility matrix -T.
 * Then the fracture bcs are imposed on -T and on the rhs.
 */
class StiffMatrixMFD: public MatrixHandlerMFD
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Construct a StiffnessMFD_Builder.
    /*!
     * @param rMesh A constant reference to the mesh
     * @param BCmap A constant reference to the boundary conditions
     * @param size  The dimension of the stiffness matrix
     */
	StiffMatrixMFD( const Rigid_Mesh & rMesh, UInt size, const BoundaryConditions & BCmap):
		MatrixHandlerMFD(rMesh, size), M_bc(BCmap), M_b(new Vector(Vector::Constant(size, 0.))){}
	//! No Copy-Constructor
    StiffMatrixMFD(const StiffMatrixMFD &) = delete;
	//! No Empty-Constructor
    StiffMatrixMFD() = delete;
	//! Destructor
    ~StiffMatrixMFD() = default;
	//@}

    //! @name Get Methods
    //@{
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
    Vector & getBCVector() const
        {return *M_b;}
    //@}

	//! @name Assemble Methods
    //@{
    //! Assemble method for the MFD stiffness matrix.
    /*!
     * Assemble the MFD stiffness matrix using the bulk and fracture builders.
     */
    void assemble();
    //@}

private:
	//! Constant reference to the boundary conditions
	const BoundaryConditions      &	M_bc;
	//! A reference to the rhs of the system
	std::unique_ptr<Vector>         M_b;  
};

} // namespace FVCode3D
#endif // __DARCYSTIFFNESS_HPP__
