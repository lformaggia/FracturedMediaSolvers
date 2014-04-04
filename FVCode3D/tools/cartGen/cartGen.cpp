//#define NDEBUG // add this macro to disable asserts
#include <cassert>
#include <iostream>
#include <chrono>
#include <vector>
#include <fstream>

#include "core/TypeDefinition.hpp"
#include "core/data.hpp"
#include "mesh/Mesh3D.hpp"
#include "mesh/Properties.hpp"
#include "mesh/cartesianGrid.hpp"
#include "import/import.hpp"
#include "export/exportVTU.hpp"
#include "utility/converter.hpp"

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
	Geometry::Mesh3D mesh;
	Geometry::PropertiesMap propMap(data.getMobility());
	std::cout << " done." << std::endl;

	std::cout << "Generate Cartesian grid..." << std::flush;
	CartesianGrid cart(mesh, propMap);
	cart.generate(true, data.getLx(), data.getLy(), data.getLz(), data.getNx(), data.getNy(), data.getNz());
	std::cout << " done." << std::endl;

	std::cout << "Compute separated cells..." << std::flush;
	mesh.updateFacetsWithCells();
	std::cout << " done." << std::endl;

	std::cout << "Compute neighboring cells..." << std::flush;
	mesh.updateCellsWithNeighbors();
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Set labels on boundary & add Fractures..." << std::flush;
	cart.extractBC();
	std::cout << " done." << std::endl;

	// create here the map "facetIdToZone"
	// Don't use "1" as zone code, it is reserved to the matrix!
	// Example:
	std::map<UInt,UInt> facetIdToZone;
	for(std::map<UInt,Geometry::Mesh3D::Facet3D>::const_iterator it = mesh.getFacetsMap().begin(); it != mesh.getFacetsMap().end(); ++it)
	{
		Geometry::Point3D centroid = it->second.getCentroid();
		if( centroid[0] > 0 && centroid[0] < 2 &&
			centroid[1] > 0 && centroid[1] < 1 &&
			centroid[2] > 0 && centroid[2] < 1)
		{
			facetIdToZone.insert( std::pair<UInt,UInt>(it->first,2));
		}
	}
	// then call "addFractures()"
	cart.addFractures(facetIdToZone);

	std::cout << "Compute facet ids of the fractures..." << std::flush;
	mesh.updateFacetsWithFractures();
	std::cout << " done." << std::endl << std::endl;

	propMap.setPropertiesOnMatrix(mesh, data.getMatrixPorosity(), data.getMatrixPermeability());
	propMap.setPropertiesOnFractures(mesh, data.getFractureAperture(), data.getFracturePorosity(), data.getFracturePermeability());

	std::cout << "Export..." << std::flush;
	saveAsSolverFormat(data.getOutputDir() + data.getOutputFile() + ".fvg", mesh, propMap);
	std::cout << " done." << std::endl;

	return 0;
}
