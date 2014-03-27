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
#include "boundaryCondition/BC.hpp"
#include "quadrature/Quadrature.hpp"
#include "import/import.hpp"
#include "export/exportVTU.hpp"
#include "utility/converter.hpp"
#include "solver/solver.hpp"
#include "problem/problem.hpp"
#include "problem/darcySteady.hpp"
#include "problem/darcyPseudoSteady.hpp"
#include "functions.hpp"

typedef Problem<EigenCholesky, CentroidQuadrature, CentroidQuadrature> Pb;
typedef DarcySteady<EigenCholesky, CentroidQuadrature, CentroidQuadrature> DarcyPb;
typedef DarcyPseudoSteady<EigenCholesky, CentroidQuadrature, CentroidQuadrature, TimeScheme::BDF2> PseudoDarcyPb;

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
	Geometry::PropertiesMap propMap(data.getMobility(), data.getCompressibility());
	std::cout << " done." << std::endl;


	std::cout << "Create Importer..." << std::flush;
	Importer * importer = 0;
	if(data.getMeshType() == Data::MeshFormatType::TPFA)
		importer = new ImporterTPFA(data.getMeshDir() + data.getMeshFile(), mesh, propMap);
	else if(data.getMeshType() == Data::MeshFormatType::forSolver)
		importer = new ImporterForSolver(data.getMeshDir() + data.getMeshFile(), mesh, propMap);
	std::cout << " done." << std::endl;


	std::cout << "Import grid file..." << std::flush;
	importer->import(data.fractureOn());
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Compute separated cells..." << std::flush;
	mesh.updateFacetsWithCells();
	std::cout << " done." << std::endl;


	std::cout << "Compute neighboring cells..." << std::flush;
	mesh.updateCellsWithNeighbors();
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Set labels on boundary" << std::flush;
	if(data.getMeshType() == Data::MeshFormatType::TPFA)
	{
		std::cout << " & add Fractures..." << std::flush;
		importer->addBCAndFractures(data.getTheta());
	}
	else if(data.getMeshType() == Data::MeshFormatType::forSolver)
	{
		std::cout << "..." << std::flush;
		importer->extractBC(data.getTheta());
	}
	std::cout << " done." << std::endl;

	std::cout << "Compute facet ids of the fractures..." << std::flush;
	mesh.updateFacetsWithFractures();
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;

	//propMap.setPropertiesOnMatrix(mesh, 0.25, 1000);
	//propMap.setPropertiesOnFractures(mesh, 0.1, 1, 100000);

	std::cout << "Export..." << std::flush;
	ExporterVTU exporter;
	exporter.exportMesh(mesh, data.getOutputDir() + data.getOutputFile() + "_mesh.vtu");
	exporter.exportFractures(mesh, data.getOutputDir() + data.getOutputFile() + "_fractures.vtu");
	exporter.exportMeshWithFractures(mesh, data.getOutputDir() + data.getOutputFile() + "_mesh_fracture.vtu");
	exporter.exportWithProperties(mesh, propMap, data.getOutputDir() + data.getOutputFile() + "_prop.vtu");
	exporter.exportWireframe(mesh, data.getOutputDir() + data.getOutputFile() + "_wire.vtu");
//	if(data.getMeshType() == Data::MeshFormatType::TPFA)
//		saveAsSolverFormat(data.getOutputDir() + data.getOutputFile() + "_new.fvg", mesh, propMap);
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Add BCs..." << std::flush;
	BoundaryConditions::BorderBC backBC	(1, Dirichlet, fZero );
	BoundaryConditions::BorderBC frontBC(2, Dirichlet, fOne );
	BoundaryConditions::BorderBC leftBC	(3, Neumann, fZero );
	BoundaryConditions::BorderBC rightBC(4, Neumann, fZero );
	BoundaryConditions::BorderBC upBC	(5, Neumann, fZero );
	BoundaryConditions::BorderBC downBC	(6, Neumann, fZero );

	std::vector<BoundaryConditions::BorderBC> borders;

	borders.push_back( backBC );
	borders.push_back( frontBC );
	borders.push_back( leftBC );
	borders.push_back( rightBC );
	borders.push_back( upBC );
	borders.push_back( downBC );

	BoundaryConditions BC(borders);
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Assemble rigid mesh..." << std::flush;
	Geometry::Rigid_Mesh myrmesh(mesh, propMap);
	std::cout << " done." << std::endl << std::endl;

	myrmesh.showMe();

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Export Fracture Junctures..." << std::flush;
	exporter.exportFractureJunctures(myrmesh, data.getOutputDir() + data.getOutputFile() + "_junc.vtu");
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Build problem..." << std::flush;
	Pb * darcy;
	if(data.getProblemType() == Data::ProblemType::steady)
		darcy = new DarcyPb(myrmesh, BC, SS);
	else if(data.getProblemType() == Data::ProblemType::pseudoSteady)
		darcy = new PseudoDarcyPb(myrmesh, BC, SS, data);
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Solve problem..." << std::flush;
	if(data.getProblemType() == Data::ProblemType::steady)
		darcy->assembleAndSolve();
	else if(data.getProblemType() == Data::ProblemType::pseudoSteady)
	{
		dynamic_cast<PseudoDarcyPb *>(darcy)->initialize();
		UInt iter=0;

		std::cout << std::endl << " ... initial solution, t = " << data.getInitialTime() << " ..." << std::endl;

		std::stringstream ss;
		ss << data.getOutputDir() + data.getOutputFile() + "_solution_";
		ss << iter;
		ss << ".vtu";
		exporter.exportSolution(myrmesh, ss.str(), dynamic_cast<PseudoDarcyPb *>(darcy)->getOldSolution());
		ss.str(std::string());
		ss << data.getOutputDir() + data.getOutputFile() + "_solution_f_";
		ss << iter;
		ss << ".vtu";
		exporter.exportSolutionOnFractures(myrmesh, ss.str(), dynamic_cast<PseudoDarcyPb *>(darcy)->getOldSolution());
		std::cout << "  done." << std::endl << std::endl;

		++iter;

		for(Real t = data.getInitialTime() + data.getTimeStep() ; t <= data.getEndTime(); t+=data.getTimeStep(), ++iter)
		{
			std::cout << " ... at t = " << t << " ..." << std::endl;
			darcy->assembleAndSolve();
			ss.str(std::string());
			std::cout << " Export Solution" << std::flush;
			ss << data.getOutputDir() + data.getOutputFile() + "_solution_";
			ss << iter;
			ss << ".vtu";
			exporter.exportSolution(myrmesh, ss.str(), darcy->getSolver().getSolution());
			ss.str(std::string());
			ss << data.getOutputDir() + data.getOutputFile() + "_solution_f_";
			ss << iter;
			ss << ".vtu";
			exporter.exportSolutionOnFractures(myrmesh, ss.str(), darcy->getSolver().getSolution());
			std::cout << "  done." << std::endl << std::endl;
		}
	}
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Export Solution..." << std::flush;
	exporter.exportSolution(myrmesh, data.getOutputDir() + data.getOutputFile() + "_solution.vtu", darcy->getSolver().getSolution());
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Export Solution on Fractures..." << std::flush;
	exporter.exportSolutionOnFractures(myrmesh, data.getOutputDir() + data.getOutputFile() + "_solution_f.vtu", darcy->getSolver().getSolution());
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Export Property..." << std::flush;
	exporter.exportWithProperties(myrmesh, data.getOutputDir() + data.getOutputFile() + "_APP.vtu", Aperture | Permeability | Porosity);
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;

	delete darcy;
	delete importer;
	
	chrono.stop();

	return 0;
}
