/*!
 *	@file Properties.hpp
 *	@brief Classes that handle the properties.
 */

#ifndef PROPERTIES_HPP_
#define PROPERTIES_HPP_

#include "core/TypeDefinition.hpp"
#include "property/PermeabilityFactory.hpp"

/* TODO
#	modified:   src/utility/converter.cpp
#	modified:   src/property/Properties.cpp
#	modified:   src/assembler/stiffness.cpp
#	modified:   src/export/exportVTU.cpp
#	modified:   src/import/import.cpp
 */

namespace FVCode3D
{

class Mesh3D;

//! Properties (aperture, porosity, permeability) of a single fracture/matrix.
/*!
 * @class Properties
 * This class contains the properties for a specific zone.
 */
class Properties
{
public:

	//! Empty constructor
	Properties():
		M_aperture(0.), M_porosity(0.), M_permeability(0.) {}

	//! Constructor
	/*!
	 * @param aperture Aperture of the fracture. If cell, should be = 1
	 * @param porosity Porosity of the matrix. If fracture, should be = 1
	 * @param permeability Permeability of the fracture/matrix
	 */
	Properties(const Real aperture, const Real porosity, const Real permeability):
		M_aperture(aperture), M_porosity(porosity), M_permeability(permeability) {}

	//! Set properties
	/*!
	 * Set properties of a single fracture/matrix.
	 * @param aperture Aperture of the fracture. If cell, should be = 1
	 * @param porosity Porosity of the matrix. If fracture, should be = 1
	 * @param permeability Permeability of the fracture/matrix
	 */
	void setProperties(const Real aperture, const Real porosity, const Real permeability)
		{ M_aperture = aperture; M_porosity = porosity;
		  M_permeability = permeability; };

	//! Aperture of the fracture. If cell, should be = 1
	Real M_aperture;
	//! Porosity of the matrix. If fracture, should be = 1
	Real M_porosity;
	//! Permeability of the fracture/matrix
	Real M_permeability;
};

//! Class that handles the properties for all the reservoir
/*!
 * @class PropertiesMap
 * This class implements a map that contains the properties of all zones
 */
class PropertiesMap
{
public:

	//! Empty constructor
	PropertiesMap():M_mobility(0.), M_compressibility(0.) {};

	//! Constructor
	/*!
	 * @param mobility mobility of the fluid
	 */
	PropertiesMap(const Real mobility, const Real compressibility = 0.):
		M_mobility(mobility), M_compressibility(compressibility) {};

	//! Get the number of zones
	/*!
	 * @return the number of zones
	 */
	UInt getNumberOfZone() const
		{ return M_properties.size(); };

	//! Get the mobility of the fluid
	/*!
	 * @return the mobility of the fluid
	 */
	Real getMobility() const
		{ return M_mobility; };

	//! Get the compressbility
	/*!
	 * @return the compressibility
	 */
	Real getCompressibility() const
		{ return M_compressibility; };

	//! Get the map of properties (const)
	/*!
	 * @return the map of the properties. First -> zone code; Second -> Properties
	 */
	const std::map<UInt, Properties> & getProperties() const
		{ return M_properties; };

	//! Get the map of properties
	/*!
	 * @return the map of the properties. First -> zone code; Second -> Properties
	 */
	std::map<UInt, Properties> & getProperties()
		{ return M_properties; };

	//! Get the properties of a selected zone (const)
	/*!
	 * @param zone the id of the zone code
	 * @return a reference of the related Properties
	 */
	const Properties & getProperties(const UInt zone) const
		{ return M_properties.at(zone); };

	//! Get the properties of a selected zone
	/*!
	 * @param zone the id of the zone code
	 * @return a reference of the related Properties
	 */
	Properties & getProperties(const UInt zone)
		{ return M_properties[zone]; };

	//! Set the same property on all the porous medium
	/*!
	 * @param mesh reference to a Mesh3D
	 * @param porosity Porosity of the matrix.
	 * @param permeability Permeability of the matrix
	 */
	void setPropertiesOnMatrix(const Mesh3D & mesh, const Real porosity, const Real permeability);

	//! Set the same property on all the fractures
	/*!
	 * @param mesh reference to a Mesh3D
	 * @param aperture Aperture of the fractures
	 * @param porosity Porosity of the fractures
	 * @param permeability Permeability of the fractures
	 */
	void setPropertiesOnFractures(const Mesh3D & mesh, const Real aperture, const Real porosity, const Real permeability);

	//! Set the properties of a selected zone
	/*!
	 * @param zone the id of the zone code
	 * @param prop a reference to a Properties
	 */
	void setZone(const UInt zone, const Properties & prop)
		{ M_properties[zone] = prop; };

	//! Set the mobility of the fluid
	/*!
	 * @param the mobility of the fluid
	 */
	void setMobility(const Real mobility)
		{ M_mobility = mobility; };

	//! Set the compressibility
	/*!
	 * @param the compressibility
	 */
	void setCompressibility(const Real compressibility)
		{ M_compressibility = compressibility; };

	//! Destructor
	~PropertiesMap() {};

private:

	//! Map that link each zone code with the property
	std::map<UInt, Properties> M_properties;

	//! Mobility of the fluid.
	Real M_mobility;

	//! Compressibility.
	Real M_compressibility;
};

} // namespace FVCode3D

#endif /* PROPERTIES_HPP_ */
