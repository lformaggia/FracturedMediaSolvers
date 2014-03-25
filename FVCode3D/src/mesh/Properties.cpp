/*!
 *	@file Properties.hpp
 *	@brief Classes that handle the properties.
 */

#include "mesh/Properties.hpp"
#include "mesh/Mesh3D.hpp"

namespace Geometry
{

void PropertiesMap::setPropertiesOnMatrix(const Mesh3D & mesh, const Real porosity, const Real permeability)
{
	std::set<UInt> zones;
	for(std::map<UInt,Mesh3D::Cell3D>::const_iterator it = mesh.getCellsMap().begin(); it != mesh.getCellsMap().end(); ++it)
		if(it->second.getZoneCode() > 0)
			zones.insert(it->second.getZoneCode());

	for(std::set<UInt>::const_iterator it = zones.begin(); it != zones.end(); ++it)
	{
		getProperties(*it).M_aperture = 1.;
		getProperties(*it).M_porosity = porosity;
		getProperties(*it).M_permeability = permeability;
	}
}

void PropertiesMap::setPropertiesOnFractures(const Mesh3D & mesh, const Real aperture, const Real porosity, const Real permeability)
{
	std::set<UInt> zones;
	for(std::map<UInt,Mesh3D::Facet3D>::const_iterator it = mesh.getFacetsMap().begin(); it != mesh.getFacetsMap().end(); ++it)
		if(it->second.getZoneCode() > 0)
			zones.insert(it->second.getZoneCode());

	for(std::set<UInt>::const_iterator it = zones.begin(); it != zones.end(); ++it)
	{
		getProperties(*it).M_aperture = aperture;
		getProperties(*it).M_porosity = porosity;
		getProperties(*it).M_permeability = permeability;
	}
}

} // namespace Geometry
