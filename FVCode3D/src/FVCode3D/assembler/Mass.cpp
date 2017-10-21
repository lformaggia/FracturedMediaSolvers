 /*!
 *  @file mass.cpp
 *  @brief Class for building a Mass-matrix for a finite volume method (definitions).
 */

#include <FVCode3D/assembler/Mass.hpp>
#include <FVCode3D/property/Properties.hpp>

namespace FVCode3D
{

void MassMatrixFV::assemble()
{
    M_matrixElements.reserve(this->M_mesh.getCellsVector().size() + 3 * this->M_mesh.getFractureFacetsIdsVector().size());

    Real porosity, volume;
    const Real compressibility = M_mesh.getPropertiesMap().getCompressibility();

    for (auto& cell_it : this->M_mesh.getCellsVector())
    {
        porosity = M_mesh.getPropertiesMap().getProperties(cell_it.getZoneCode()).M_porosity;
        M_matrixElements.emplace_back(M_offsetRow + cell_it.getId(), M_offsetCol + cell_it.getId(), compressibility * porosity * cell_it.getVolume()); // Triplet
    }

    for (auto& facet_it : this->M_mesh.getFractureFacetsIdsVector())
    {
        // set mass term of the fracture
        volume = M_mesh.getPropertiesMap().getProperties(facet_it.getZoneCode()).M_aperture * facet_it.getFacet().area();
        porosity = M_mesh.getPropertiesMap().getProperties(facet_it.getZoneCode()).M_porosity;

        M_matrixElements.emplace_back(M_offsetRow + facet_it.getIdAsCell(), M_offsetCol + facet_it.getIdAsCell(), compressibility * porosity * volume); // Triplet
    }

//    if( this->M_mesh.getCellsVector().size() > 0 )
//    {
//        for (auto& facet_it : this->M_mesh.getFractureFacetsIdsVector())
//        {
//            // set mass term of the fracture
//            volume = M_mesh.getPropertiesMap().getProperties(facet_it.getZoneCode()).M_aperture * facet_it.getFacet().area();
//
//            // remove from the first neighbor the portion due to the fracture
//            porosity = M_mesh.getPropertiesMap().getProperties( this->M_mesh.getCellsVector()[facet_it.getSeparatedCellsIds()[0]].getZoneCode() ).M_porosity;
//            M_matrixElements.emplace_back(M_offsetRow + facet_it.getSeparatedCellsIds()[0], M_offsetCol + facet_it.getSeparatedCellsIds()[0], -compressibility * porosity * volume / 2); // Triplet
//
//            // remove from the second neighbor the portion due to the fracture
//            porosity = M_mesh.getPropertiesMap().getProperties( this->M_mesh.getCellsVector()[facet_it.getSeparatedCellsIds()[1]].getZoneCode() ).M_porosity;
//            M_matrixElements.emplace_back(M_offsetRow + facet_it.getSeparatedCellsIds()[1], M_offsetCol + facet_it.getSeparatedCellsIds()[1], -compressibility * porosity * volume / 2); // Triplet
//        }
//    }
}

} // namespace FVCode3D
