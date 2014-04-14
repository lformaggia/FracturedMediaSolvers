 /*!
 *	@file mass.cpp
 *	@brief Class for building a Mass-matrix for a finite volume method (definitions).
 */

#include "assembler/mass.hpp"
#include "mesh/Properties.hpp"

namespace Darcy
{

void MassMatrix::assemble()
{
	std::vector<Triplet> Matrix_elements;
	Real porosity, volume;
	const Real compressibility = M_properties.getCompressibility();

	for (auto cell_it : this->M_mesh.getCellsVector())
	{
		porosity = M_properties.getProperties(cell_it.getZoneCode()).M_porosity;
		Matrix_elements.emplace_back(Triplet (cell_it.getId(), cell_it.getId(), compressibility * porosity * cell_it.getVolume()));
	}

	for (auto facet_it : this->M_mesh.getFractureFacetsIdsVector())
	{
		// set mass term of the fracture
		volume = M_properties.getProperties(facet_it.getZoneCode()).M_aperture * facet_it.getFacet().area();
		porosity = M_properties.getProperties(facet_it.getZoneCode()).M_porosity;

		Matrix_elements.emplace_back(Triplet (facet_it.getIdasCell(), facet_it.getIdasCell(), compressibility * porosity * volume));

		// remove from the first neighbor the portion due to the fracture
		porosity = M_properties.getProperties( this->M_mesh.getCellsVector()[facet_it.getSeparated()[0]].getZoneCode() ).M_porosity;
		Matrix_elements.emplace_back(Triplet (facet_it.getSeparated()[0], facet_it.getSeparated()[0], -compressibility * porosity * volume / 2));

		// remove from the second neighbor the portion due to the fracture
		porosity = M_properties.getProperties( this->M_mesh.getCellsVector()[facet_it.getSeparated()[1]].getZoneCode() ).M_porosity;
		Matrix_elements.emplace_back(Triplet (facet_it.getSeparated()[1], facet_it.getSeparated()[1], -compressibility * porosity * volume / 2));
	}
	this->M_Matrix->setFromTriplets(Matrix_elements.begin(), Matrix_elements.end());
}

} // namespace Darcy
