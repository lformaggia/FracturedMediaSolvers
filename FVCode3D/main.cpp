//#define NDEBUG // add this macro to disable asserts
#include <cassert>
#include <iostream>
#include <chrono>
#include <vector>
#include <fstream>
#include <functional>

#include "core/TypeDefinition.hpp"
#include "mesh/Rigid_Mesh.hpp"
#include "mesh/Properties.hpp"
#include "assembler/stiffness.hpp"
#include "assembler/mass.hpp"
#include "boundaryCondition/BC.hpp"
#include "quadrature/Quadrature.hpp"
#include "import/import.hpp"
#include "export/export.hpp"
#include "utility/converter.hpp"

typedef	std::function<Real (Geometry::Point3D point)> Func;

int main()
{
	//defining the path to save results
	std::string res_folder("./results/");
	std::string grid_file("data/grid2.grid");
	std::string output_file("grid2");

	Chrono chrono;
	chrono.start();

	std::cout << "Define Mesh and Properties..." << std::flush;
	Geometry::Mesh3D mesh;
	Geometry::PropertiesMap propMap(2.27);
	std::cout << " done." << std::endl;


	std::cout << "Create Importer..." << std::flush;
	ImporterTPFAStandard importer(grid_file, mesh, propMap);
	//ImporterTPFAWithBC import("data/test.grid", mesh, propMap);
	std::cout << " done." << std::endl;


	std::cout << "Import grid file..." << std::flush;
	importer.import();
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Compute separated cells..." << std::flush;
	mesh.updateFacetsWithCells();
	std::cout << " done." << std::endl;


	std::cout << "Compute neighboring cells..." << std::flush;
	mesh.updateCellsWithNeighbors();
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Set BC on boundary & add Fractures..." << std::flush;
	readTPFAStandardAsTPFAWithBC(mesh, propMap);
	std::cout << " done." << std::endl;

	std::cout << "Compute facet ids of the fractures..." << std::flush;
	mesh.updateFacetsWithFractures();
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Export..." << std::flush;
	ExporterVTU exporter;
	exporter.exportMesh(mesh, res_folder + output_file + "_mesh.vtu");
	exporter.exportFractures(mesh, res_folder + output_file + "_fractures.vtu");
	exporter.exportMeshWithFractures(mesh, res_folder + output_file + "_mesh_fracture.vtu");
	convertFromTPFAStandardToTPFAWithBC(res_folder + output_file + "_bc.grid", mesh, propMap);
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Add BC..." << std::flush;
	Darcy::BoundaryConditions::BorderBC upBC (5, Darcy::Neumann, [](Geometry::Point3D) { return 0.; } );
	Darcy::BoundaryConditions::BorderBC downBC (6, Darcy::Neumann, [](Geometry::Point3D) { return 0.; } );
	Darcy::BoundaryConditions::BorderBC backBC (1, Darcy::Neumann, [](Geometry::Point3D) { return 0.; } );
	Darcy::BoundaryConditions::BorderBC frontBC (2, Darcy::Neumann, [](Geometry::Point3D) { return 0.; } );
	Darcy::BoundaryConditions::BorderBC leftBC (3, Darcy::Neumann, [](Geometry::Point3D) { return -1.; } );
	Darcy::BoundaryConditions::BorderBC rightBC (4, Darcy::Neumann, [](Geometry::Point3D) { return 1.; } );

	std::vector<Darcy::BoundaryConditions::BorderBC> borders;

	borders.push_back( upBC );
	borders.push_back( downBC );
	borders.push_back( backBC );
	borders.push_back( frontBC );
	borders.push_back( leftBC );
	borders.push_back( rightBC );

	Darcy::BoundaryConditions BC(borders);
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Assembling rigid mesh..." << std::flush;
	Geometry::Rigid_Mesh myrmesh(mesh, propMap);
	std::cout << " done." << std::endl << std::endl;

	myrmesh.showMe();

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Assembling stiffness matrix..." << std::flush;
	Darcy::StiffMatrix S(myrmesh, BC);
	S.assemble();
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


//	std::cout << "Assembling mass matrix..." << std::flush;
//	Darcy::MassMatrix M(myrmesh);
//	M.assemble();
//	std::cout << " done." << std::endl << std::endl;
//
//	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Assembling source vector..." << std::flush;
	//auto SorceDomain = [](Geometry::Point3D p){return (p.x()*p.x() + p.y()*p.y() + p.z()*p.z()) < 1;};
	//auto SinkDomain = [](Geometry::Point3D p){return (p.x()*p.x() + p.y()*p.y()) < 1;};
	auto Source = [](Geometry::Point3D p){return 100*( ( (p.x()-3.15389e+07)*(p.x()-3.15389e+07) +
												     (p.y()-1.68937e+07)*(p.y()-1.68937e+07) +
												     (p.z()+14048.1    )*(p.z()+14048.1    ) )
												     <=5e4
												 ); };
	auto Sink = [](Geometry::Point3D p){return -100*( ( (p.x()-3.14986e+07)*(p.x()-3.14986e+07) +
		     	 	 	 	 	 	 	 	 	    (p.y()-1.68939e+07)*(p.y()-1.68939e+07) +
		     	 	 	 	 	 	 	 	 	    (p.z()+14041.5    )*(p.z()+14041.5    ) )
		     	 	 	 	 	 	 	 	 	    <=5e4
		 	 	 	 	 	 	 	 	 	   ); };
//	auto Source = [](Geometry::Point3D p){return 100*( ( (p.x()-3.15185e7)*(p.x()-3.15185e7) +
//												     (p.y()-1.6913e7)*(p.y()-1.6913e7) +
//												     (p.z()+1.80e4    )*(p.z()+1.80e4    ) )
//												     <=5e4
//												 ); };
//	auto Sink = [](Geometry::Point3D p){return -100*( ( (p.x()-3.15185e7)*(p.x()-3.15185e7) +
//		     	 	 	 	 	 	 	 	 	    (p.y()-1.6874e7)*(p.y()-1.6874e7) +
//		     	 	 	 	 	 	 	 	 	    (p.z()+1.635e4    )*(p.z()+1.635e4    ) )
//		     	 	 	 	 	 	 	 	 	    <=5e4
//		 	 	 	 	 	 	 	 	 	   ); };
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
	b = S.getBCVector();
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Solve problem..." << std::flush;

	Eigen::SimplicialCholesky <SpMat> chol(A);
	Eigen::VectorXd x = chol.solve(b);
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Export Solution..." << std::flush;
	exporter.exportSolution(myrmesh, res_folder + output_file + "_solution.vtu", x);
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Export Solution on Fractures..." << std::flush;
	exporter.exportSolutionOnFractures(myrmesh, res_folder + output_file + "_solution_f.vtu", x);
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


	std::cout << "Export Property..." << std::flush;
	exporter.exportWithProperties(myrmesh, res_folder + output_file + "_perm.vtu", Permeability);
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Export Property..." << std::flush;
	exporter.exportWithProperties(myrmesh, res_folder + output_file + "_poro.vtu", Porosity);
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;

	chrono.stop();

	return 0;
}
