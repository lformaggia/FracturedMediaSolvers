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
#include "assembler/MatrixHandler.hpp"
#include "boundaryCondition/BC.hpp"

class Rigid_mesh;
class PropertiesMap;

namespace Darcy
{

//! Class for assembling a stiffness matrix
/*!
    @class StiffMatrix
    This class constructs the stiffness-matrix for the Darcy problem.
    The adopted technique is a two point finite volume method.
    The fractures are considered as cells and take part to discretization.
*/
class StiffMatrix: public MatrixHandler
{

    //! Typedef for std::pair<UInt,UInt>
    /*!
        @typedef Fracture_Juncture
        This type definition permits to treat std::pair<UInt,UInt> as a Fracture_Juncture.
    */
    typedef std::pair<UInt,UInt> Fracture_Juncture;
    //! Typedef for Geometry::Rigid_Mesh::Facet_ID
    /*!
        @typedef Facet_ID
        This type definition permits to treat Geometry::Rigid_Mesh::Facet_ID as a Facet_ID.
    */
    typedef Geometry::Rigid_Mesh::Facet_ID Facet_ID;
    //! Typedef for Geometry::Rigid_Mesh::Edge_ID
    /*!
        @typedef Edge_ID
        This type definition permits to treat Geometry::Rigid_Mesh::Edge_ID as a Edge_ID.
    */
    typedef Geometry::Rigid_Mesh::Edge_ID Edge_ID;

public:
    //! @name Constructor & Destructor
    //@{

    //! Construct a stiffness-Matrix, given a Geometry::Rigid_Mesh and the boundary conditions
    /*!
        @param rigid_mesh A Geometry::Rigid_Mesh used to build the matrix
        @param BC Boundary conditions given in the container Darcy::BoundaryConditions
    */
    StiffMatrix(const Geometry::Rigid_Mesh & rigid_mesh, const BoundaryConditions & Bc):
        MatrixHandler(rigid_mesh), _b (new Vector(Vector::Constant( this->M_size, 0.))),
        M_properties(rigid_mesh.getPropertiesMap()), m_Bc(Bc) {}
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
        {return *_b;}
    //@}

    //! @name Methods
    //@{
    //! Assemble method
    /*!
     * @return Assemble the stiffness matrix
     */
    void assemble();
    //@}

public:

	//! @name Alpha Methods
	//@{

	//! Border center
	/*!
	 * @param fj the Id of the juncture of two Fracture_Facet in 3D
	 * @return The center of the juncture between two Fracure_Facet
	 */
	Generic_Point border_center(Fracture_Juncture fj) const;

	//! It is called by the method assemble() and it computes the coefficient alpha
	/*!
	 * @param cellId the Id of a Cell
	 * @param facet A pointer to a Geometry::Rigid_mesh::Facet_ID
	 * @return The computed coefficient alpha
	 */
	Real Findalpha (const UInt & cellId, const Facet_ID * facet) const;

	//! It is called by the method assemble() and it computes the coefficient alpha
	/*!
	 * @param facetId the Id of a Facet
	 * @param edge A pointer to a Geometry::Rigid_mesh::Edge_ID
	 * @return The computed coefficient alpha
	 */
	Real Findalpha (const UInt & facetId, const Edge_ID * edge) const;

	//! It is called by the method assemble() and it computes the coefficient alpha in the case of Dirichlet BC
	/*!
	 * @param cellId the Id of a Cell
	 * @param facet A pointer to a Geometry::Rigid_mesh::Facet_ID
	 * @return The computed coefficient alpha
	 */
	Real FindDirichletalpha (const UInt & cellId, const Facet_ID * facet) const;

	//! It is called by the method assemble() and it computes the coefficient alpha in the case of Dirichlet BC
	/*!
	 * @param facetId the Id of a Facet
	 * @param edge A pointer to a Geometry::Rigid_mesh::Edge_ID
	 * @return The computed coefficient alpha
	 */
	Real FindDirichletalpha (const UInt & facetId, const Edge_ID * edge) const;

	//! It is called by the method assemble() and it computes the coefficient alpha in the case of a fracture in 3D
	/*!
	 * @param fj is a Fracture_Juncture
	 * @param n_Id The Id of the Fracture_Facet
	 * @return The computed coefficient alpha
	 */
	Real Findfracturesalpha (const std::pair<UInt,UInt> fj, const UInt n_Id) const;
	//@}

protected:    
    
    //! @name Assemble Methods
    //@{
    
    //! Assemble the porous medium block
    /*!
     * @return Assemble the porous media block
     */
    void assemblePorousMatrix( std::vector<Triplet>& Matrix_elements ) const;
    
    //! Assemble the BCs for the porous medium
    /*!
     * @return Assemble the BCs for the porous medium
     */
    void assemblePorousMatrixBC( std::vector<Triplet>& Matrix_elements ) const;

    //! Assemble the fractures block
    /*!
     * @return Assemble the fractures block
     */
    void assembleFracture( std::vector<Triplet>& Matrix_elements ) const;
    
    //! Assemble the BCs for the fractures
    /*!
     * @return Assemble the BCs for the fractures
     */
    void assembleFractureBC( std::vector<Triplet>& Matrix_elements ) const;
    //@}
    
protected:
    //! Unique pointer to the vector that contains the effects of BCs on the RHS
    std::unique_ptr<Vector> _b;
    //! A reference to a Geometry::PropertiesMap
    const Geometry::PropertiesMap & M_properties;
    //! The container of the BCs
    const BoundaryConditions & m_Bc;
};

} // namespace Darcy

#endif
