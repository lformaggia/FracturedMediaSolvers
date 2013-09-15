
#include<iostream>
#include<chrono>
#include<vector>

#include "helperfunction.hpp"
#include "GetPot"
#include "./src_Turconi/Fracture2D.hpp"
#include "./src_Turconi/Mesh2D.hpp"
#include "./src_Turconi/Triangulation2D.hpp"
#include "./DP_src/Rigid_Mesh.hpp"
#include "./DP_Darcy/DarcyStiffness.hpp"
#include "./DP_Darcy/DarcyMass.hpp"
#include "./DP_Darcy/MatrixHandler.hpp"
#include "./DP_Darcy/DarcyBC.hpp"
#include "./DP_Darcy/DarcyQuadrature.hpp"
#include "./DP_Darcy/DarcyQuadratureRules.hpp"


int main(int argc, char** argv)
{

	using namespace std::chrono;
	typedef	std::function<double (Geometry::Point2D point)> Func;
	typedef Geometry::Point2D Point;
	typedef Geometry::Dimension<2> R2;
	auto start = high_resolution_clock::now();

	//read parameters from input
	double L_max = 0.07;
	Func u, derxu, deryu;
	readParameters(argc,argv, L_max, u, derxu, deryu);
	auto Source = [&u](Geometry::Point2D point){return (3.*u(point));};

	// define domain border
	double pi=4.*std::atan(1.0);
	std::vector<Geometry::Point2D> borderVertexes;

	borderVertexes.push_back( Geometry::Point2D(pi,pi) );
	borderVertexes.push_back( Geometry::Point2D(pi,0) );
	borderVertexes.push_back( Geometry::Point2D(0,0) );
	borderVertexes.push_back( Geometry::Point2D(0,pi) );

	// define DomainBorder
	std::vector<Geometry::Domain<R2>::DomainBorder> domainborders;

	Geometry::Domain<R2>::DomainBorder border0 (Geometry::Segment2D (borderVertexes[0], borderVertexes[1]));
	Geometry::Domain<R2>::DomainBorder border1 (Geometry::Segment2D (borderVertexes[1], borderVertexes[2]));
	Geometry::Domain<R2>::DomainBorder border2 (Geometry::Segment2D (borderVertexes[2], borderVertexes[3]));
	Geometry::Domain<R2>::DomainBorder border3 (Geometry::Segment2D (borderVertexes[3], borderVertexes[0]));

	//build Domain
	domainborders.push_back ( border0 );
	domainborders.push_back ( border1 );
	domainborders.push_back ( border2 );
	domainborders.push_back ( border3 );

	Geometry::Domain<Geometry::Dimension<2> > domain(domainborders);


	//define Boundary Conditions
	Darcy::BoundaryConditions<R2>::BorderBC BC1 (0, Darcy::Dirichlet, [&u](Point p){return u(p);});
	Darcy::BoundaryConditions<R2>::BorderBC BC2 (1, Darcy::Neumann, [&deryu, &pi](Point p){return -deryu(Point(p.x(),0));});			
	Darcy::BoundaryConditions<R2>::BorderBC BC3 (2, Darcy::Dirichlet, [&u](Point p){return u(p);});		
	Darcy::BoundaryConditions<R2>::BorderBC BC4 (3, Darcy::Neumann, [&deryu, &pi](Point p){return deryu(Point(p.x(),pi));});

	//build BoundaryConditions
	std::vector<Darcy::BoundaryConditions<R2>::BorderBC> borders;

	borders.push_back ( BC1 );
	borders.push_back ( BC3 );
	borders.push_back ( BC4 );
	borders.push_back ( BC2 );

	Darcy::BoundaryConditions<R2> BC(borders, domain);

	//defining the path to save vtk
	std::string vtk_folder("./Darcy_vtk/");

	// default triangulation
	Geometry::Triangulation2D Triang(borderVertexes);
	Triang.buildTriangulation(20.705, L_max);

	//assembling mesh2D
	std::cout << std::endl;
	auto step0 = high_resolution_clock::now();
	std::cout << duration_cast<seconds>(step0-start).count() << "seconds are passed" << std::endl;
	std::cout << "Assembling mesh2D" << std::endl;

	Geometry::Mesh2D mesh;
	mesh.addTriangulation(Triang);

	//assembling Rigid_Mesh
	std::cout << std::endl;
	auto step1 = high_resolution_clock::now();
	std::cout << duration_cast<seconds>(step1-start).count() << " seconds are passed" << std::endl;
	std::cout  << "Assembling rigid mesh" << std::endl;

	Geometry::Rigid_Mesh<Geometry::Dimension<2> > myrmesh(mesh,domain);

	//some information about the mesh	
	std::cout << std::endl;
	std::cout <<"numero di nodi: "<<myrmesh.getNodesVector().size()<<" ,mesh2D: "<<mesh.getNodesVector().size()<<std::endl;
	std::cout<<"numero di lati: "<<myrmesh.getFacetsVector().size()<<" ,mesh2D: "<<mesh.getEdgesSet().size()<<std::endl;
	std::cout<<"numero di celle: "<<myrmesh.getCellsVector().size()<<" ,mesh2D: "<<mesh.getCellsMap().size()<<std::endl;
	std::cout << std::endl;

	//some information on the Facets size
	std::cout << "Max Facet lenght: " << myrmesh.getMaxFacetSize () << std::endl;		
	std::cout << "Min Facet lenght: " << myrmesh.getMinFacetSize () << std::endl;		
	std::cout << "Average Facet lenght: " << myrmesh.getAveFacetSize () << std::endl;		
	std::cout << std::endl;


	//create stiffness matrix A
	auto step2 = high_resolution_clock::now();
	std::cout << duration_cast<seconds>(step2-start).count() << " seconds are passed" << " : " <<std::endl;
	std::cout <<"Assembling stiffness matrix"<<std::endl;
	std::cout << std::endl;

	Darcy::StiffMatrix<Geometry::Dimension<2> > S(myrmesh, BC, [](Geometry::Point2D point){return 1.;},1.);
	S.assemble();

	//create mass matrix B
	auto step3 = high_resolution_clock::now();
	std::cout << duration_cast<seconds>(step3-start).count() << " seconds are passed" << std::endl;
	std::cout <<"Assembling mass matrix"<<std::endl;
	std::cout << std::endl;

	Darcy::MassMatrix<Geometry::Dimension<2> > M(myrmesh);
	M.assemble();

	//define problem Matrix M
	Darcy::SpMat A = S.getMatrix() + M.getMatrix();

	//compute the known term
	Darcy::Quadrature<Geometry::Dimension<2> > quadrature2(myrmesh, Darcy::Triangle2D());

	Eigen::VectorXd f(S.getSize());
	f = quadrature2.CellIntegrate (Source);

	Eigen::VectorXd g(S.getSize());
	g = S.getBCVector() + f;	

	//discretize the exact solution
	Eigen::VectorXd uh (S.getSize());
	uh = quadrature2.CellIntegrate (u);
	Eigen::SimplicialCholesky <Darcy::SpMat> inve(M.getMatrix());
	Eigen::VectorXd Uh = inve.solve(uh);

	//solve the problem
	auto step4 = high_resolution_clock::now();
	std::cout << duration_cast<seconds>(step4-start).count() << " seconds are passed" << std::endl;
	std::cout <<"Solving problem"<<std::endl;
	std::cout << std::endl;

	Eigen::SimplicialCholesky <Darcy::SpMat> chol(A);
	Eigen::VectorXd x = chol.solve(g);

	//saving mesh and solution in vtk
	myrmesh.exportVtk(vtk_folder+"TriangDefault.vtk");
	myrmesh.appendSolutionToVtk(x, vtk_folder+"TriangDefault.vtk");

	//compute the L2 norm of the difference exact_solution-computed_solution
	Eigen::VectorXd solDiff (S.getSize());
	solDiff = x-Uh;
	double L2_Norm = quadrature2.L2Norm(solDiff);
	std::cout << "NORMA L2 DELL'ERRORE: " << L2_Norm << std::endl;
	std::cout << std::endl;

	//save the errorrs
	std::fstream filestr;
	std::string errorsfile = "L2Errors";
	filestr.open (errorsfile.c_str(), std::ios_base::out | std::ios_base::app	);
	if (filestr.is_open())
	{
		std::cout << std::endl << " File: " << errorsfile << ", successfully opened";
	}
	else
	{
		std::cerr << std::endl << " *** Error: file not opened *** " << std::endl << std::endl;
	}
	std::cout << std::endl << " Appending L2 errore to the file... " << std::endl << std::endl;

	filestr << myrmesh.getMaxFacetSize() << "\t" << myrmesh.getAveFacetSize() << "\t" << L2_Norm;
	filestr << std::endl;

	//clear Mesh2d
	mesh.clear();

	//final time
	auto step5 = high_resolution_clock::now();
	std::cout << duration_cast<seconds>(step5-start).count() << " seconds is the total time" << std::endl;
	std::cout << std::endl;

	return 0;
}
