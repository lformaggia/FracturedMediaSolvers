//#define NDEBUG // add this macro to disable asserts
#include <cassert>
#include <iostream>
#include <chrono>
#include <vector>
#include <fstream>

#include "core/TypeDefinition.hpp"
#include "core/data.hpp"
#include "mesh/Rigid_Mesh.hpp"
#include "mesh/Properties.hpp"
#include "assembler/stiffness.hpp"
#include "assembler/mass.hpp"
#include "boundaryCondition/BC.hpp"
#include "quadrature/Quadrature.hpp"
#include "import/import.hpp"
#include "export/export.hpp"
#include "utility/converter.hpp"
#include "solver/solver.hpp"
#include "functions.hpp"

int main(int argc, char * argv[])
{
	GetPot command_line(argc,argv);
	const std::string dataFileName = command_line.follow("data.txt", 2, "-f", "--file");

	Chrono chrono;
	chrono.start();

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
//	ImporterTPFA importer(data.getMeshDir() + data.getMeshFile(), mesh, propMap);
	ImporterForSolver importer(data.getMeshDir() + data.getMeshFile(), mesh, propMap);
	std::cout << " done." << std::endl;


	std::cout << "Import grid file..." << std::flush;
	importer.import(data.fractureOn());
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Compute separated cells..." << std::flush;
	mesh.updateFacetsWithCells();
	std::cout << " done." << std::endl;


	std::cout << "Compute neighboring cells..." << std::flush;
	mesh.updateCellsWithNeighbors();
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Set labels on boundary & add Fractures..." << std::flush;
	importer.addBCAndFractures(data.getTheta());
	std::cout << " done." << std::endl;

	std::cout << "Compute facet ids of the fractures..." << std::flush;
	mesh.updateFacetsWithFractures();
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Export..." << std::flush;
	ExporterVTU exporter;
//	exporter.exportMesh(mesh, data.getOutputDir() + data.getOutputFile() + "_mesh.vtu");
//	exporter.exportFractures(mesh, data.getOutputDir() + data.getOutputFile() + "_fractures.vtu");
//	exporter.exportMeshWithFractures(mesh, data.getOutputDir() + data.getOutputFile() + "_mesh_fracture.vtu");
//	saveAsSolverFormat(data.getOutputDir() + data.getOutputFile() + "_bc.grid", mesh, propMap);
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Add BCs..." << std::flush;
	Darcy::BoundaryConditions::BorderBC backBC	(1, Darcy::Neumann, fZero );
	Darcy::BoundaryConditions::BorderBC frontBC	(2, Darcy::Neumann, fZero );
	Darcy::BoundaryConditions::BorderBC leftBC	(3, Darcy::Dirichlet, fOne );
	Darcy::BoundaryConditions::BorderBC rightBC	(4, Darcy::Dirichlet, fMinusOne );
	Darcy::BoundaryConditions::BorderBC upBC	(5, Darcy::Neumann, fZero );
	Darcy::BoundaryConditions::BorderBC downBC	(6, Darcy::Neumann, fZero );

	std::vector<Darcy::BoundaryConditions::BorderBC> borders;

	borders.push_back( backBC );
	borders.push_back( frontBC );
	borders.push_back( leftBC );
	borders.push_back( rightBC );
	borders.push_back( upBC );
	borders.push_back( downBC );

	Darcy::BoundaryConditions BC(borders);
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Assemble rigid mesh..." << std::flush;
	Geometry::Rigid_Mesh myrmesh(mesh, propMap);
	std::cout << " done." << std::endl << std::endl;

	myrmesh.showMe();

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Export Fracture Junctures..." << std::flush;
//	exporter.exportFractureJunctures(myrmesh, data.getOutputDir() + data.getOutputFile() + "_junc.vtu");
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Assemble stiffness matrix..." << std::flush;
	Darcy::StiffMatrix S(myrmesh, BC);
	S.assemble();
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Assemble source/sink term..." << std::flush;

	Darcy::Quadrature quad(myrmesh, Darcy::CentroidQuadrature(), Darcy::CentroidQuadrature());
	Eigen::VectorXd f1(S.getSize());
	Eigen::VectorXd f2(S.getSize());
	Eigen::VectorXd f(S.getSize());
	f1 = quad.CellIntegrate(Source);
	f2 = quad.CellIntegrate(Sink);
	f = f1 + f2;
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Define problem..." << std::flush;
	SpMat A = S.getMatrix() ;
	Eigen::VectorXd b;
	b = S.getBCVector() + f;
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Solve problem..." << std::flush;

	EigenCholesky solver(A,b);
	solver.solve();
	const Vector & x = solver.getSolution();
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Export Solution..." << std::flush;
	exporter.exportSolution(myrmesh, data.getOutputDir() + data.getOutputFile() + "_solution.vtu", x);
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Export Solution on Fractures..." << std::flush;
	exporter.exportSolutionOnFractures(myrmesh, data.getOutputDir() + data.getOutputFile() + "_solution_f.vtu", x);
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Export Property..." << std::flush;
	exporter.exportWithProperties(myrmesh, data.getOutputDir() + data.getOutputFile() + "_perm.vtu", Permeability);
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Export Property..." << std::flush;
	exporter.exportWithProperties(myrmesh, data.getOutputDir() + data.getOutputFile() + "_poro.vtu", Porosity);
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;

	chrono.stop();

	return 0;
}
