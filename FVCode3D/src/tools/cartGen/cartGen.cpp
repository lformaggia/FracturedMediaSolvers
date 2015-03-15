//#define NDEBUG // add this macro to disable asserts
#include <cassert>
#include <iostream>
#include <chrono>
#include <vector>
#include <fstream>

#include <FVCode3D/core/TypeDefinition.hpp>
#include <FVCode3D/core/Data.hpp>
#include <FVCode3D/mesh/Mesh3D.hpp>
#include <FVCode3D/property/Properties.hpp>
#include <FVCode3D/property/Permeability.hpp>
#include <FVCode3D/mesh/CartesianGrid.hpp>
#include <FVCode3D/import/Import.hpp>
#include <FVCode3D/export/ExportVTU.hpp>
#include <FVCode3D/export/ExportCP.hpp>
#include <FVCode3D/utility/Converter.hpp>

using namespace FVCode3D;

int main(int argc, char * argv[])
{
    GetPot command_line(argc,argv);
    const std::string dataFileName = command_line.follow("data.txt", 2, "-f", "--file");

    std::cout << "Read Data..." << std::flush;
    Data data(dataFileName);
    data.verbose(true);
    std::cout << " done." << std::endl;

    std::cout << std::endl;
    data.showMe();
    std::cout << std::endl;

    std::cout << "Define Mesh and Properties..." << std::flush;
    Mesh3D mesh;
    PropertiesMap propMap(data.getMobility());
    std::cout << " done." << std::endl;

    std::cout << "Generate Cartesian grid..." << std::flush;
    CartesianGrid cart(mesh, propMap);
    cart.generate(true, data.getLx(), data.getLy(), data.getLz(), data.getNx(), data.getNy(), data.getNz(),
                        data.getSx(), data.getSy(), data.getSz(),
                        data.noiseOn(), data.getMeanNormalDistribution(), data.getStDevNormalDistribution()); // 5.0e-2  7.5e-3 , 3.75e-3 , 1.875e-3 , 9.375e-4
    std::cout << " done." << std::endl << std::endl;

    std::cout << "# of cells: " << mesh.getCellsMap().size() << std::endl << std::endl;

    std::cout << "Compute separated cells..." << std::flush;
    mesh.updateFacetsWithCells();
    std::cout << " done." << std::endl;

    std::cout << "Compute neighboring cells..." << std::flush;
    mesh.updateCellsWithNeighbors();
    std::cout << " done." << std::endl << std::endl;

    std::cout << "Set labels on boundary & add Fractures..." << std::flush;
    cart.extractBC();
    std::cout << " done." << std::endl << std::endl;

    // create here the map "facetIdToZone"
    // Don't use "1" as zone code, it is reserved to the matrix!
    // Example:
    std::map<UInt,UInt> facetIdToZone;
    /*
    for(std::map<UInt,Mesh3D::Facet3D>::const_iterator it = mesh.getFacetsMap().begin(); it != mesh.getFacetsMap().end(); ++it)
    {
        Point3D centroid = it->second.getCentroid();
        // along x
        for(Real thresh = 0.1; thresh<2 ; thresh+=0.1)
        {
            if(centroid[0] > 1. + 0.05 || centroid[0] < 1. - 0.05)
            {
                if( centroid[0] > thresh-0.01 && centroid[0] < thresh + 0.01 &&
                    centroid[1] > 0.06 && centroid[1] < 0.94 &&
                    centroid[2] > 0.06 && centroid[2] < 0.94)
                {
                    facetIdToZone.insert( std::pair<UInt,UInt>(it->first,2));
                }
            }
        }
        // along y
        for(Real thresh = 0.1; thresh<1; thresh+=0.1)
        {
            if( centroid[0] > 0.06 && centroid[0] < 1.94 &&
                centroid[1] > thresh-0.01 && centroid[1] < thresh + 0.01 &&
                centroid[2] > 0.06 && centroid[2] < 0.94)
            {
                facetIdToZone.insert( std::pair<UInt,UInt>(it->first,2));
            }
        }
        // along z
        for(Real thresh = 0.1; thresh<1; thresh+=0.1)
        {
            if( centroid[0] > 0.06 && centroid[0] < 1.94 &&
                centroid[1] > 0.06 && centroid[1] < 0.94 &&
                centroid[2] > thresh-0.01 && centroid[2] < thresh + 0.01)
            {
                facetIdToZone.insert( std::pair<UInt,UInt>(it->first,2));
            }
        }
    }*/
    /*
    for(std::map<UInt,Mesh3D::Facet3D>::const_iterator it = mesh.getFacetsMap().begin(); it != mesh.getFacetsMap().end(); ++it)
    {
        Point3D centroid = it->second.getCentroid();
        if(
            (
                (centroid[1] > 0. - 0.01 && centroid[1] < 0. + 0.01)
                &&
                (centroid[0] > -0.4 - 0.01 && centroid[0] < 0.8 + 0.01)
            )
            ||
            (
                (centroid[0] > 0.4 - 0.01 && centroid[0] < 0.8 + 0.01)
                &&
                (centroid[2] > -0.25 - 0.01 && centroid[2] < -0.25 + 0.01)
                &&
                (centroid[1] > -0.4 - 0.01 && centroid[1] < 0.7 + 0.01)
            )
            ||
            (
                (centroid[2] > 0.75 - 0.01 && centroid[2] < 0.75 + 0.01)
                &&
                (centroid[1] > -0.25 - 0.01 && centroid[1] < 0.25 + 0.01)
                &&
                (centroid[0] > -0.5 - 0.01 && centroid[0] < 0.75 + 0.01)
            )
            ||
            (
                (centroid[2] > -0.75 - 0.01 && centroid[2] < -0.75 + 0.01)
                &&
                (centroid[1] > -0.75 - 0.01 && centroid[1] < 0.25 + 0.01)
                &&
                (centroid[0] > -0.5 - 0.01 && centroid[0] < -0.25 + 0.01)
            )
            ||
            (
                (centroid[2] > 0. - 0.01 && centroid[2] < 0. + 0.01)
                &&
                (centroid[1] > -1. - 0.01 && centroid[1] < -0.25 + 0.01)
                &&
                (centroid[0] > -0.75 - 0.01 && centroid[0] < 0.75 + 0.01)
            )

          )
        {
            facetIdToZone.insert( std::pair<UInt,UInt>(it->first,2));
        }
    }
    */
    // then call "addFractures()"
    cart.addFractures(facetIdToZone);

    std::cout << "# of fracture facets: " << facetIdToZone.size() << std::endl << std::endl;

    std::cout << "Compute facet ids of the fractures..." << std::flush;
    mesh.updateFacetsWithFractures();
    std::cout << " done." << std::endl << std::endl;

    std::cout << "Set uniform properties..." << std::flush;
    std::shared_ptr<PermeabilityBase> matrixPerm( new PermeabilityScalar );
    std::shared_ptr<PermeabilityBase> fracturesPerm( new PermeabilityScalar );
    matrixPerm->setPermeability( data.getMatrixPermeability(), 0 );
    fracturesPerm->setPermeability( data.getFracturePermeability(), 0 );
    propMap.setPropertiesOnMatrix(mesh, data.getMatrixPorosity(), matrixPerm);
    propMap.setPropertiesOnFractures(mesh, data.getFractureAperture(), data.getFracturePorosity(), fracturesPerm);
    std::cout << " done." << std::endl << std::endl;

    std::cout << "Export..." << std::flush;
    saveAsSolverFormat(data.getOutputDir() + data.getOutputFile() + ".fvg", mesh, propMap);

    ExporterCP exporter;
    exporter.exportMesh(mesh, data.getOutputDir() + data.getOutputFile() + ".GRDECL");
    std::cout << " done." << std::endl;

    return 0;
}
