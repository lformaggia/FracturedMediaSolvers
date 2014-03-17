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
#include "import/import.hpp"
#include "export/export.hpp"
#include "utility/converter.hpp"

int main(int argc, char * argv[])
{
	GetPot command_line(argc,argv);
	const std::string dataFileName = command_line.follow("data.txt", 2, "-f", "--file");

	std::cout << "Read Data..." << std::flush;
	Data data(dataFileName);
	std::cout << " done." << std::endl;

	std::cout << std::endl;
	data.showMe();
	std::cout << std::endl;

	std::cout << "Define Mesh and Properties..." << std::flush;
	Geometry::Mesh3D mesh;
	Geometry::PropertiesMap propMap(data.getMobility());
	std::cout << " done." << std::endl;

	std::cout << "Create Importer..." << std::flush;
	ImporterTPFA importer(data.getMeshDir() + data.getMeshFile(), mesh, propMap);
	std::cout << " done." << std::endl;

	std::cout << "Import grid file..." << std::flush;
	importer.import(data.fractureOn());
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Compute separated cells..." << std::flush;
	mesh.updateFacetsWithCells();
	std::cout << " done." << std::endl;

	std::cout << "Compute neighboring cells..." << std::flush;
	mesh.updateCellsWithNeighbors();
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Set labels on boundary & add Fractures..." << std::flush;
	importer.addBCAndFractures(data.getTheta());
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Compute facet ids of the fractures..." << std::flush;
	mesh.updateFacetsWithFractures();
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Save as solver format..." << std::flush;
	saveAsSolverFormat(data.getOutputDir() + data.getOutputFile() + ".grid", mesh, propMap);
	std::cout << " done." << std::endl << std::endl;

	return 0;
}