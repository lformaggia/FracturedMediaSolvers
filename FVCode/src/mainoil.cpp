
#include<iostream>
#include<chrono>
#include<vector>

#include "./src_Turconi/Fracture2D.hpp"
#include "./src_Turconi/Mesh2D.hpp"
#include "./src_Turconi/Triangulation2D.hpp"
#include "./DP_src/Rigid_Mesh.hpp"
#include "./DP_Darcy/DarcyStiffness.hpp"
#include "./DP_Darcy/DarcyMass.hpp"
#include "./DP_Darcy/DarcyBC.hpp"
#include "./DP_Darcy/DarcyQuadrature.hpp"
#include "./DP_Darcy/DarcyQuadratureRules.hpp"

//#include "./my_src/Dimension.hpp"


int main()
{
	using namespace std::chrono;
	system_clock::time_point tp = system_clock::now();
	system_clock::duration dtn = tp.time_since_epoch();


	// define domain border
	std::vector<Geometry::Point2D> borderVertexes;
	
	borderVertexes.push_back( Geometry::Point2D(12,10) );
	borderVertexes.push_back( Geometry::Point2D(12,0) );
	borderVertexes.push_back( Geometry::Point2D(0,0) );
	borderVertexes.push_back( Geometry::Point2D(0,10) );

	// define DomainBorder
	std::vector<Geometry::Domain<Geometry::Dimension<2> >::DomainBorder> domainborders;

	Geometry::Domain<Geometry::Dimension<2> >::DomainBorder border0 (Geometry::Segment2D (borderVertexes[0], borderVertexes[1]));
	Geometry::Domain<Geometry::Dimension<2> >::DomainBorder border1 (Geometry::Segment2D (borderVertexes[1], Geometry::Point2D(1,0)));
	Geometry::Domain<Geometry::Dimension<2> >::DomainBorder border12 (Geometry::Segment2D (Geometry::Point2D(1,0), borderVertexes[2]));
	Geometry::Domain<Geometry::Dimension<2> >::DomainBorder border2 (Geometry::Segment2D (borderVertexes[2], borderVertexes[3]));
	Geometry::Domain<Geometry::Dimension<2> >::DomainBorder border23 (Geometry::Segment2D (borderVertexes[3], Geometry::Point2D(11,10)));
	Geometry::Domain<Geometry::Dimension<2> >::DomainBorder border3 (Geometry::Segment2D (Geometry::Point2D(11,10), borderVertexes[0]));

	//build Domain
	domainborders.push_back ( border0 );
	domainborders.push_back ( border1 );
	domainborders.push_back ( border12 );
	domainborders.push_back ( border2 );
	domainborders.push_back ( border23 );
	domainborders.push_back ( border3 );

	Geometry::Domain<Geometry::Dimension<2> > domain(domainborders);

	//define Boundary Conditions
	Darcy::BoundaryConditions<Geometry::Dimension<2> >::BorderBC BC1 (0, Darcy::Neumann, [](Geometry::Point2D ){return 0.;});
	Darcy::BoundaryConditions<Geometry::Dimension<2> >::BorderBC BC12 (1, Darcy::Dirichlet, [](Geometry::Point2D ){return 10.;});
	Darcy::BoundaryConditions<Geometry::Dimension<2> >::BorderBC BC2(2, Darcy::Neumann, [](Geometry::Point2D ){return 10.;});			
	Darcy::BoundaryConditions<Geometry::Dimension<2> >::BorderBC BC3 (3, Darcy::Neumann, [](Geometry::Point2D ){return 0.;});		
	Darcy::BoundaryConditions<Geometry::Dimension<2> >::BorderBC BC34 (4, Darcy::Dirichlet, [](Geometry::Point2D ){return 0.;});
	Darcy::BoundaryConditions<Geometry::Dimension<2> >::BorderBC BC4 (5, Darcy::Neumann, [](Geometry::Point2D ){return -10.;});

	//build BoundaryConditions
	std::vector<Darcy::BoundaryConditions<Geometry::Dimension<2> >::BorderBC> borders;

	borders.push_back ( BC1 );
	borders.push_back ( BC12 );
	borders.push_back ( BC2 );
	borders.push_back ( BC3 );
	borders.push_back ( BC34 );
	borders.push_back ( BC4 );

	Darcy::BoundaryConditions<Geometry::Dimension<2> > BC(borders, domain);


	// define fracture
	Geometry::Fracture2D f1( Geometry::Point2D(2,8),
						  Geometry::Point2D(9,8) );
	f1.setPermeability(1.e+3);
	f1.setAperture(0.1);
	Geometry::Fracture2D f2( Geometry::Point2D(2,8),
						  Geometry::Point2D(2,4) );
	f2.setPermeability(1.e+3);
	f2.setAperture(0.1);
	Geometry::Fracture2D f3( Geometry::Point2D(1.5,4),
						  Geometry::Point2D(4.5,7) );
	f3.setPermeability(1.e+3);
	f3.setAperture(0.1);
	Geometry::Fracture2D f4( Geometry::Point2D(3,7),
						  Geometry::Point2D(5,5) );
	f4.setPermeability(1.e+3);
	f4.setAperture(0.1);
	Geometry::Fracture2D f5( Geometry::Point2D(5,3),
						  Geometry::Point2D(6,3) );
	f5.setPermeability(1.e+3);
	f5.setAperture(0.1);
	Geometry::Fracture2D f6( Geometry::Point2D(4,2),
						  Geometry::Point2D(10,6) );
	f6.setPermeability(1.e+3);
	f6.setAperture(0.1);
	Geometry::Fracture2D f7( Geometry::Point2D(9,4),
						  Geometry::Point2D(8,5) );
	f7.setPermeability(1.e+3);
	f7.setAperture(0.1);


	// build fracture network
	Geometry::FractureNetwork2D fn;
	fn.push_back(f1);
	fn.push_back(f2);
	fn.push_back(f3);
	fn.push_back(f4);
	fn.push_back(f5);
	fn.push_back(f6);
	fn.push_back(f7);

	//defining the path to save vtk
	std::string vtk_folder("./Darcy_vtk/");

	// adaptive triangulation
	std::cout << std::endl;
	std::cout <<"Building adaptive triangulation"<<std::endl;
	system_clock::time_point tp2 = system_clock::now();
	system_clock::duration dtn2 = tp2.time_since_epoch();
	std::cout << "Passed seconds: " << (dtn2.count() - dtn.count())* system_clock::period::num / system_clock::period::den;
	std::cout << std::endl;
	std::cout << std::endl;

	Geometry::Triangulation2D Triang(borderVertexes);

	Triang.addFractureNetwork(fn);
	Triang.buildAdaptiveTriangulation();
	Triang.buildAdaptiveTriangulation(20.705,0.1,5.0,3.0,1.0);
	Triang.findApproxFractureNetwork();

	//assembling mesh2D
	std::cout << std::endl;
	std::cout <<"Assembling mesh2D"<<std::endl;
	tp2 = system_clock::now();
	dtn2 = tp2.time_since_epoch();
	std::cout << "Passed seconds: " << (dtn2.count() - dtn.count())* system_clock::period::num / system_clock::period::den;
	std::cout << std::endl;

	Geometry::Mesh2D mesh;
	mesh.addTriangulation(Triang);

	//assembling Rigid_Mesh
	std::cout << std::endl;
	std::cout <<"Assembling rigid mesh"<<std::endl;
	tp2 = system_clock::now();
	dtn2 = tp2.time_since_epoch();
	std::cout << "Passed seconds: " << (dtn2.count() - dtn.count())* system_clock::period::num / system_clock::period::den;
	std::cout << std::endl;
	std::cout << std::endl;

	Geometry::Rigid_Mesh<Geometry::Dimension<2> > myrmesh(mesh,domain);

	//some information about the mesh	
	std::cout<<"numero di nodi: "<<myrmesh.getNodesVector().size()<<" ,mesh2D: "<<mesh.getNodesVector().size()<<std::endl;
	std::cout<<"numero di lati: "<<myrmesh.getFacetsVector().size()<<" ,mesh2D: "<<mesh.getEdgesSet().size()<<std::endl;
	std::cout<<"numero di celle: "<<myrmesh.getCellsVector().size()<<" ,mesh2D: "<<mesh.getCellsMap().size()<<std::endl;
	std::cout << std::endl;

	myrmesh.showMe();
	std::cout << std::endl;
	myrmesh.getFacetsVector()[0].showMe();
	std::cout << std::endl;
	myrmesh.getCellsVector().rbegin()->showMe();
	std::cout << std::endl;

	//create stiffness matrix A
	std::cout <<"Assembling stiffness matrix"<<std::endl;
	tp2 = system_clock::now();
	dtn2 = tp2.time_since_epoch();
	std::cout << "Passed seconds: " << (dtn2.count() - dtn.count())* system_clock::period::num / system_clock::period::den;
	std::cout << std::endl;
	std::cout << std::endl;

	Darcy::StiffMatrix<Geometry::Dimension<2> > A(myrmesh, BC, [](Geometry::Point2D point){return 2.54;}, 2.27);
	A.assemble();


	//create mass matrix B
	std::cout <<"Assembling mass matrix"<<std::endl;
	tp2 = system_clock::now();
	dtn2 = tp2.time_since_epoch();
	std::cout << "Passed seconds: " << (dtn2.count() - dtn.count())* system_clock::period::num / system_clock::period::den;
	std::cout << std::endl;
	std::cout << std::endl;

	Darcy::MassMatrix<Geometry::Dimension<2> > B(myrmesh);
	B.assemble();

	
	//define problem Matrix M
	Darcy::SpMat M = A.getMatrix() ;

	Eigen::VectorXd g = A.getBCVector();

	//solve the problem
	std::cout <<"solving problem"<<std::endl;
	tp2 = system_clock::now();
	dtn2 = tp2.time_since_epoch();
	std::cout << "Passed seconds: " << (dtn2.count() - dtn.count())* system_clock::period::num / system_clock::period::den;
	std::cout << std::endl;

	Eigen::SimplicialCholesky <Darcy::SpMat> chol(M);
	Eigen::VectorXd x = chol.solve(g);

	//saving mesh and solution in vtk
	myrmesh.exportVtk(vtk_folder+"TriangDefault.vtk");
	myrmesh.appendSolutionToVtk(x, vtk_folder+"TriangDefault.vtk");

	//clear Mesh2d
	mesh.clear();

	return 0;
}
