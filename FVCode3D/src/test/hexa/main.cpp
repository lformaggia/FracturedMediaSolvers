// Add this macro to disable asserts
//#define NDEBUG

// Add this macro the enable the exporter
//#define FVCODE3D_EXPORT

#include <cassert>
#include <iostream>
#include <chrono>
#include <vector>
#include <fstream>
#include <cmath>

#include <FVCode3D/FVCode3D.hpp>

using namespace FVCode3D;

typedef Problem<CentroidQuadrature, CentroidQuadrature> Pb;
typedef DarcySteady<CentroidQuadrature, CentroidQuadrature> DarcyPb;

int main(int argc, char * argv[])
{
    Chrono chrono;
    chrono.start();

    GetPot command_line(argc,argv);
    const std::string dataFileName = command_line.follow("data.txt", 2, "-f", "--file");

    DataPtr_Type dataPtr(new Data(dataFileName));

    Mesh3D mesh;
    PropertiesMap propMap(dataPtr->getMobility(), dataPtr->getCompressibility());

    Importer * importer = 0;
    importer = new ImporterForSolver(dataPtr->getMeshDir() + dataPtr->getMeshFile(), mesh, propMap);

    importer->import(dataPtr->fractureOn());

    mesh.updateFacetsWithCells();
    mesh.updateCellsWithNeighbors();

    importer->extractBC(dataPtr->getTheta());

    mesh.updateFacetsWithFractures();

    std::shared_ptr<PermeabilityBase> matrixPerm( new PermeabilityScalar );
    std::shared_ptr<PermeabilityBase> fracturesPerm( new PermeabilityScalar );
    const Real kf = 1.e-3;
    matrixPerm->setPermeability( 1., 0 );
    fracturesPerm->setPermeability( kf, 0 );
    const Real aperture = 1.e-2;
    const Real matrixPoro = 0.25;
    const Real fracturesPoro = 1.;
    propMap.setPropertiesOnMatrix(mesh, matrixPoro, matrixPerm);
    propMap.setPropertiesOnFractures(mesh, aperture, fracturesPoro, fracturesPerm);

    ExporterVTU exporter;

	Func uex1 = [kf,aperture](Point3D p)
		{ return  kf*sin(p.x())*cosh(p.y()) + (1.-kf)/2.*cosh(aperture/2.)*sin(p.x())*cos(p.z()); };
	
	Func uex2 = [kf](Point3D p)
		{ return  -kf*cos(p.x())*sinh(p.y()); };
	
	Func uex3 = [kf,aperture](Point3D p)
		{ return  (1.-kf)/2.*cosh(aperture/2.)*cos(p.x())*sin(p.z()); };

	Func pex = [kf,aperture](Point3D p)
		{ return  kf*cos(p.x())*cosh(p.y()) + (1.-kf)/2.*cosh(aperture/2.)*cos(p.x())*cos(p.z()); };

	Func pfex = [kf,aperture](Point3D p)
		{ return  kf*cos(p.x()) + (1.-kf)/2.*cosh(aperture/2.)*cos(p.x())*cos(p.z()); };

	Func SS = [kf,aperture](Point3D p)
		{ return  p.y() == 0 ? 	
			pow(kf,2)*cos(p.x()) + kf*(1.-kf)*cosh(aperture/2.)*cos(p.x())*cos(p.z())
			:
			(1.-kf)*cosh(aperture/2.)*cos(p.x())*cos(p.z()); };
			
	Func fZero = [](Point3D){ return 0.; };

    BoundaryConditions::BorderBC leftBC (BorderLabel::Left, Dirichlet, pex );
    BoundaryConditions::BorderBC rightBC(BorderLabel::Right, Dirichlet, pex );
    BoundaryConditions::BorderBC backBC (BorderLabel::Back, Dirichlet, pex );
    BoundaryConditions::BorderBC frontBC(BorderLabel::Front, Dirichlet, pex );
    BoundaryConditions::BorderBC upBC   (BorderLabel::Top, Dirichlet, pex );
    BoundaryConditions::BorderBC downBC (BorderLabel::Bottom, Dirichlet, pex );

    std::vector<BoundaryConditions::BorderBC> borders;

    borders.push_back( backBC );
    borders.push_back( frontBC );
    borders.push_back( leftBC );
    borders.push_back( rightBC );
    borders.push_back( upBC );
    borders.push_back( downBC );

    BoundaryConditions BC(borders);

    Rigid_Mesh myrmesh(mesh, propMap);
    myrmesh.BuildEdges();
    myrmesh.BuildOrientationFacets();

    myrmesh.showMe();
    std::cout << "hMin: " << myrmesh.getMinEdgeSize() << std::endl;
    std::cout << "hMax: " << myrmesh.getMaxEdgeSize() << std::endl;
    std::cout << "hAve: " << myrmesh.getAveEdgeSize() << std::endl;
	std::cout << std::endl;
    
    Pb * darcy(nullptr);
    darcy = new DarcyPb(dataPtr->getSolverType(), myrmesh, BC, SS, dataPtr);

    darcy->assemble();

	if(dynamic_cast<IterativeSolver*>(darcy->getSolverPtr()))
    {
		dynamic_cast<IterativeSolver*>(darcy->getSolverPtr())->setMaxIter(dataPtr->getIterativeSolverMaxIter());
		dynamic_cast<IterativeSolver*>(darcy->getSolverPtr())->setTolerance(dataPtr->getIterativeSolverTolerance());
		dynamic_cast<IterativeSolver*>(darcy->getSolverPtr())->set_precon(dataPtr->getpreconType());
	}

    Real Tsolving1 = chrono.partial();
    darcy->solve();
    if(dynamic_cast<IterativeSolver*>(darcy->getSolverPtr()))
    {
		dynamic_cast<IterativeSolver*>(darcy->getSolverPtr())->print();
        std::cout<<std::endl;
	}
	std::cout<<"Time to solve the system: "<< chrono.partial()-Tsolving1 <<"s."<<std::endl<<std::endl;

	UInt numFacetsTot   = myrmesh.getFacetsVector().size() + myrmesh.getFractureFacetsIdsVector().size();
	UInt numCellsTot    = myrmesh.getCellsVector().size() + myrmesh.getFractureFacetsIdsVector().size();
	std::cout << "Export Solution..." << std::flush;
	exporter.exportSolution( myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_solution.vtu", 
		darcy->getSolver().getSolution().tail(numCellsTot) );
	std::cout << "Export Solution on Fractures..." << std::flush;
	exporter.exportSolutionOnFractures( myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_solution_f.vtu", 
		darcy->getSolver().getSolution().tail(numCellsTot) );
	std::cout << " done." << std::endl << std::endl;

	UInt numCells       = myrmesh.getCellsVector().size();
	UInt numFracs       = myrmesh.getFractureFacetsIdsVector().size();

	Vector u   = darcy->getSolver().getSolution().head(numFacetsTot);
	Vector p   = darcy->getSolver().getSolution().tail(numCellsTot);

	Vector ut1 = Vector::Zero(numFacetsTot);
	Vector ut2 = Vector::Zero(numFacetsTot);
	Vector ut3 = Vector::Zero(numFacetsTot);

	Vector utn  = Vector::Zero(numFacetsTot);
	Vector pt   = evaluate(myrmesh,pex);

	for(auto & facet_it : myrmesh.getFacetsVector())
	{
		ut1[facet_it.getId()] = uex1(facet_it.getCentroid());
		ut2[facet_it.getId()] = uex2(facet_it.getCentroid());
		ut3[facet_it.getId()] = uex3(facet_it.getCentroid());
		Point3D normal = facet_it.getUnsignedNormal();
		utn[facet_it.getId()] = ut1[facet_it.getId()]*normal.x() + ut2[facet_it.getId()]*normal.y() + ut3[facet_it.getId()]*normal.z();
	}

	Vector p_diff     = pt - p;
	Vector u_diff     = utn - u;

	Quadrature q(myrmesh);
	const Real err_p    = q.L2NormMatrix(p_diff)/q.L2NormMatrix(pt);
	const Real err_pf   = q.L2NormFractures(p_diff)/q.L2NormFractures(pt);
	
	global_InnerProduct gIP(myrmesh,dFacet,dFacet);
	gIP.reserve_space();
	gIP.assemble();
	const Real norm_vel = sqrt(utn.transpose()*gIP.getMatrix()*utn);
	const Real err_u    = sqrt(u_diff.transpose()*gIP.getMatrix()*u_diff)/norm_vel;
	

	std::cout << std::setprecision(15)<<"Pressure error = "<<err_p<<std::endl<<std::endl;
	std::cout << std::setprecision(15)<<"Pressure in fracture error = "<<err_pf<<std::endl<<std::endl;
	std::cout << std::setprecision(15)<<"Velocity error = "<<err_u<<std::endl<<std::endl;


    delete darcy;
    delete importer;

    chrono.stop();
    return 0;
}

