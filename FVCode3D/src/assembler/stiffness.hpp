/*!
 *	@file stiffness.hpp
 *	@brief This class build a Stiffness-matrix of the Darcy problem.
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

	typedef std::pair<UInt,UInt> Fracture_Juncture;
	typedef Geometry::Rigid_Mesh::Facet_ID Facet_ID;

public:
	//! @name Constructor & Destructor
	//@{

	//! Construct a stiffness-Matrix, given a Geometry::Rigid_Mesh and the boundary conditions
	/*!
		@param rigid_mesh A Geometry::Rigid_Mesh used to build the matrix
		@param BC Boundary conditions given in the container Darcy::BoundaryConditions
	*/
	StiffMatrix(const Geometry::Rigid_Mesh & rigid_mesh, BoundaryConditions & Bc):
		MatrixHandler(rigid_mesh), _b (new Vector(this->M_size)),
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
	 * @return Assemble the Mass matrix
	 */
	void assemble();
	//@}

protected:

	//! @name Protected Methods
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
	Real Findalpha (const UInt & cellId, Facet_ID * const facet) const;

	//! It is called by the method assemble() and it computes the coefficient alpha in the case of Dirichlet BC
	/*!
	 * @param cellId the Id of a Cell
	 * @param facet A pointer to a Geometry::Rigid_mesh::Facet_ID
	 * @return The computed coefficient alpha
	 */
	Real FindDirichletalpha (const UInt & cellId, Facet_ID * const facet) const;

	//! It is called by the method assemble() and it computes the coefficient alpha in the case of a fracture in 3D
	/*!
	 * @param fj is a Fracture_Juncture
	 * @param n_Id The Id of the Fracture_Facet
	 * @return The computed coefficient alpha
	 */
	Real Findfracturesalpha (const std::pair<UInt,UInt> fj, const UInt n_Id) const;
	//@}
	
protected:
	//! Unique pointer to the vector that contains the effects of BCs on the RHS
	std::unique_ptr<Vector> _b;
	//! A reference to a Geometry::PropertiesMap
	const Geometry::PropertiesMap & M_properties;
	//! The container of the BCs
	BoundaryConditions & m_Bc;
};

} // namespace Darcy

#endif