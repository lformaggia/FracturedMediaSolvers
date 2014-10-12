/*!
 *  @file stiffness.hpp
 *  @brief This class build a Stiffness-matrix of the Darcy problem.
 */

#ifndef __DARCYSTIFFNESS_HPP__
#define __DARCYSTIFFNESS_HPP__

#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
#include <FVCode3D/assembler/MatrixHandler.hpp>
#include <FVCode3D/boundaryCondition/BC.hpp>

namespace FVCode3D
{

class Rigid_mesh;
class PropertiesMap;

//! Class for assembling a stiffness matrix
/*!
 * @class StiffMatrix
 * This class constructs the stiffness-matrix for the Darcy problem.
 * The adopted technique is a two point finite volume method.
 * The fractures are considered as cells and take part to discretization.
 */
class StiffMatrix: public MatrixHandler
{

    //! Typedef for std::pair<UInt,UInt>
    /*!
        @typedef Fracture_Juncture
        This type definition permits to treat std::pair<UInt,UInt> as a Fracture_Juncture.
    */
    typedef std::pair<UInt,UInt> Fracture_Juncture;
    //! Typedef for Rigid_Mesh::Facet_ID
    /*!
        @typedef Facet_ID
        This type definition permits to treat Rigid_Mesh::Facet_ID as a Facet_ID.
    */
    typedef Rigid_Mesh::Facet_ID Facet_ID;
    //! Typedef for Rigid_Mesh::Edge_ID
    /*!
        @typedef Edge_ID
        This type definition permits to treat Rigid_Mesh::Edge_ID as a Edge_ID.
    */
    typedef Rigid_Mesh::Edge_ID Edge_ID;

public:
    //! @name Constructor & Destructor
    //@{

    //! Construct a stiffness-Matrix, given a Rigid_Mesh and the boundary conditions
    /*!
        @param rigid_mesh A Rigid_Mesh used to build the matrix
        @param BC Boundary conditions given in the container BoundaryConditions
    */
    StiffMatrix(const Rigid_Mesh & rigid_mesh, const BoundaryConditions & BC):
        MatrixHandler(rigid_mesh), M_b (new Vector(Vector::Constant( this->M_size, 0.))),
        M_bc(BC) {}
    //! No Copy-Constructor
    StiffMatrix(const StiffMatrix &) = delete;
    //! No Empty-Constructor
    StiffMatrix() = delete;
    //! Destructor
    ~StiffMatrix() = default;
    //@}

    //! @name Get Methods
    //@{
    //! Get BC vector (const)
    /*!
     * @return A reference to a constant vector that represents the part of the right hand side due to the boundary conditions.
     */
    const Vector & getBCVector() const
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

} // namespace FVCode3D
#endif // __DARCYSTIFFNESS_HPP__
