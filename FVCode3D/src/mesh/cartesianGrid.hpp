/*!
 * @file cartesianGrid.hpp
 * @brief Class that generate a hexahedral structured (Cartesian) grid.
 */

#ifndef CARTESIANGRID_HPP_
#define CARTESIANGRID_HPP_

#include "core/TypeDefinition.hpp"

namespace FVCode3D
{

class Mesh3D;
class PropertiesMap;

//! Class used to generate hexahedaral structured meshes.
/*!
 * @class CartesianGrid
 * This class that allows to generate Cartesian grids.
 */
class CartesianGrid
{
public:

	//! Constructor
	/*!
	 * @param mesh reference to a Geometry::Mesh3D
	 * @param properties reference to a Geometry::PropertiesMap
	 */
	CartesianGrid(Geometry::Mesh3D & mesh, Geometry::PropertiesMap & properties):
		M_mesh(mesh), M_properties(properties) {}

	//! Generate a Cartesian mesh
	/*!
	 * @param Lx domain dimension along x axis
	 * @param Ly domain dimension along y axis
	 * @param Lz domain dimension along z axis
	 * @param Nx number of elements along x axis
	 * @param Ny number of elements along y axis
	 * @param Nz number of elements along z axis
	 */
	virtual void generate(bool fracturesOn = true, const Real Lx = 2., const Real Ly = 1., const Real Lz = 1., const UInt Nx = 40, const UInt Ny = 20, const UInt Nz = 20);

	//! Generate the BC ids
	/*!
	 * Add the BCs ids to the boundary facets.
	 * The id is set by considering the maximum component and the sign of the normal of a boundary facet.
	 * @param theta rotation angle along z-axis. It is used only to compute the BC ids. Default = 0
	 * @pre generate()
	 */
	virtual void extractBC(const Real theta = 0.);

	//! Add the fractures network
	/*!
	 * Add the fracture network from a list of facet ids
	 * @param facetIdToZone map that associates to a facet id the zone code
	 * @pre generate()
	 */
	virtual void addFractures(const std::map<UInt,UInt> & facetIdToZone);

	//! Add the lacking data
	/*!
	 * Add the BCs and the fracture network
	 * @param facetIdToZone map that associates to a facet id the zone code
	 * @param theta rotation angle along z-axis. It is used only to compute the BC ids. Default = 0
	 * @pre generate()
	 */
	virtual void addBCAndFractures(const std::map<UInt,UInt> & facetIdToZone, const Real theta = 0.);

	//! Get mesh
	/*!
	 * @return reference to the mesh
	 */
	const Geometry::Mesh3D & getMesh() const
		{ return M_mesh; }

	//! Get properties
	/*!
	 * @return reference to the properties
	 */
	const Geometry::PropertiesMap & getProperties() const
		{ return M_properties; }

	//! Destructor
	virtual ~CartesianGrid() = default;

protected:

	//! Reference to a Geometry::Mesh3D
	Geometry::Mesh3D & M_mesh;
	//! Reference to a Geometry::PropertiesMap
	Geometry::PropertiesMap & M_properties;

private:

	//! No default constructor
	CartesianGrid();

	//! No copy-constructor
	CartesianGrid(const CartesianGrid &);

	//! No assignment operator
	CartesianGrid & operator=(const CartesianGrid &);

};

} // FVCode3D
#endif /* CARTESIANGRID_HPP_ */
