/*!
 *	@file mass.hpp
 *	@brief This class build a Mass-matrix for a finite volume method.
 */ 

#ifndef __DARCYMASS_HPP__
#define __DARCYMASS_HPP__

#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
#include "assembler/MatrixHandler.hpp"
	
class PropertiesMap;

namespace FVCode3D
{

//! Class for assembling a mass matrix
/*!
	@class MassMatrix
	This class constructs the mass-matrix.
	The adopted technique is the one of the finite volume method:
	it hence represents the volume of the cell. The fractures are considered as cells.
*/
class MassMatrix: public MatrixHandler
{

	typedef std::pair<UInt,UInt> Fracture_Juncture;

public:
	//! @name Constructor & Destructor
	//@{

	//! Construct a Mass-Matrix, given a Geometry::Rigid_Mesh.
	/*!
		@param rigid_mesh A Geometry::Rigid_Mesh used to build the matrix
	*/
	MassMatrix(const Geometry::Rigid_Mesh & rigid_mesh):
		MatrixHandler(rigid_mesh, D_Cell), M_properties(rigid_mesh.getPropertiesMap()) {}
	//! No Copy-Constructor
	MassMatrix(const MassMatrix&) = delete;
	//! No Empty-Constructor
	MassMatrix() = delete;
	//! Destructor
	~MassMatrix() = default;
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
	//! A reference to a Geometry::PropertiesMap
	const Geometry::PropertiesMap & M_properties;
};

} //namespace FVCode3D

#endif
