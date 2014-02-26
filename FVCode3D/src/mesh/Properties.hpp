/*!
 *	@file Properties.hpp
 *	@brief Classes that handle the properties.
 */

#ifndef PROPERTIES_HPP_
#define PROPERTIES_HPP_

#include "core/TypeDefinition.hpp"

namespace Geometry{

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
	PropertiesMap():M_mobility(0.) {};

	//! Constructor
	/*!
	 * @param mobility mobility of the fluid
	 */
	PropertiesMap(const Real mobility):
		M_mobility(mobility) {};

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

	//! Get the map of properties (const)
	/*!
	 * @return the map of the properties. First -> zone code; Second -> Geometry::Properties
	 */
	const std::map<UInt, Properties> & getProperties() const
		{ return M_properties; };

	//! Get the map of properties
	/*!
	 * @return the map of the properties. First -> zone code; Second -> Geometry::Properties
	 */
	std::map<UInt, Properties> & getProperties()
		{ return M_properties; };

	//! Get the properties of a selected zone (const)
	/*!
	 * @param zone the id of the zone code
	 * @return a reference of the related Geometry::Properties
	 */
	const Properties & getProperties(const UInt zone) const
		{ return M_properties.at(zone); };

	//! Get the properties of a selected zone
	/*!
	 * @param zone the id of the zone code
	 * @return a reference of the related Geometry::Properties
	 */
	Properties & getProperties(const UInt zone)
		{ return M_properties[zone]; };

	//! Set the properties of a selected zone
	/*!
	 * @param zone the id of the zone code
	 * @param prop a reference to a Geometry::Properties
	 */
	void setZone(const UInt zone, const Properties & prop)
		{ M_properties[zone] = prop; };

	//! Set the mobility of the fluid
	/*!
	 * @param the mobility of the fluid
	 */
	void setMobility(const Real mobility)
		{ M_mobility =  mobility; };

	//! Destructor
	~PropertiesMap() {};

private:

	//! Map that link each zone code with the property
	std::map<UInt, Properties> M_properties;

	//! Mobility of the fluid.
	Real M_mobility;
};

} //namespace Geometry

#endif /* PROPERTIES_HPP_ */
