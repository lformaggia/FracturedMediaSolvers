// Add this macro to disable asserts
//#define NDEBUG

// Add this macro the enable the exporter
//#define FVCODE3D_EXPORT

#include <cassert>
#include <iostream>
#include <chrono>
#include <vector>
#include <fstream>

#include <unsupported/Eigen/SparseExtra>
#include <FVCode3D/core/TypeDefinition.hpp>
#include <FVCode3D/core/Data.hpp>
#include <FVCode3D/mesh/RigidMesh.hpp>
#include <FVCode3D/property/Properties.hpp>
#include <FVCode3D/property/Permeability.hpp>
#include <FVCode3D/preconditioner/preconditioner.hpp>
#include <FVCode3D/mesh/CartesianGrid.hpp>
#include <FVCode3D/boundaryCondition/BC.hpp>
#include <FVCode3D/quadrature/Quadrature.hpp>
#include <FVCode3D/import/Import.hpp>
#include <FVCode3D/export/ExportVTU.hpp>
#include <FVCode3D/utility/Converter.hpp>
#include <FVCode3D/problem/Problem.hpp>
#include <FVCode3D/problem/DarcySteady.hpp>
#include <FVCode3D/problem/DarcyPseudoSteady.hpp>
#include <FVCode3D/assembler/FixPressureDofs.hpp>
#include <FVCode3D/utility/Evaluate.hpp>
#include <FVCode3D/utility/readFunctionsFromGetpot.hpp>
#include <FVCode3D/utility/dumpMatrices.hpp>

using namespace FVCode3D;


typedef Problem<CentroidQuadrature, CentroidQuadrature> Pb;
typedef DarcySteady<CentroidQuadrature, CentroidQuadrature> DarcyPb;
typedef DarcyPseudoSteady<CentroidQuadrature, CentroidQuadrature, TimeScheme::BDF2> PseudoDarcyPb;

int main(int argc, char * argv[])
{
    GetPot command_line(argc,argv);
    const std::string dataFileName = command_line.follow("data.txt", 2, "-f", "--file");

    Chrono chrono;
    chrono.start();

    std::cout << "Read Data from " << dataFileName << "..." << std::flush;
    //DataPtr_Type dataPtr(new Data(dataFileName));
    // use global variable
    dataPtr->load(dataFileName);
    std::cout << " done." << std::endl;

    // Now for boundary conditionst and source term
    FVCode3D::ReadFunctionsFromGetpot functions(dataFileName);

    std::cout << std::endl;
    dataPtr->showMe();
    std::cout << std::endl;


    std::cout << "Define Mesh and Properties..." << std::flush;
    Mesh3D mesh;
    PropertiesMap propMap(dataPtr->getMobility(), dataPtr->getCompressibility());
    std::cout << " done." << std::endl;


    std::cout << "Create Importer..." << std::flush;
    Importer * importer = 0;
    if(dataPtr->getMeshType() == Data::MeshFormatType::TPFA)
        importer = new ImporterTPFA(dataPtr->getMeshDir() + dataPtr->getMeshFile(), mesh, propMap);
    else if(dataPtr->getMeshType() == Data::MeshFormatType::forSolver)
        importer = new ImporterForSolver(dataPtr->getMeshDir() + dataPtr->getMeshFile(), mesh, propMap);
    std::cout << " done." << std::endl;


    std::cout << "Import grid file..." << std::flush;
    importer->import(dataPtr->fractureOn());
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
    if(dataPtr->getMeshType() == Data::MeshFormatType::TPFA)
    {
        std::cout << " & add Fractures..." << std::flush;
        importer->addBCAndFractures(dataPtr->getTheta());
    }
    else if(dataPtr->getMeshType() == Data::MeshFormatType::forSolver)
    {
        std::cout << "..." << std::flush;
        importer->extractBC(dataPtr->getTheta());
    }
    std::cout << " done." << std::endl;

    std::cout << "Compute facet ids of the fractures..." << std::flush;
    mesh.updateFacetsWithFractures();
    std::cout << " done." << std::endl << std::endl;

    std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;

    std::cout << "Set properties..." << std::flush;
    /* Old version

    std::shared_ptr<PermeabilityBase> matrixPerm( new PermeabilityDiagonal );
    std::shared_ptr<PermeabilityBase> fracturesPerm( new PermeabilityDiagonal );
//	std::shared_ptr<PermeabilityBase> fracturesPerm( new PermeabilityScalar );
    matrixPerm->setPermeability( 1., 0 );
    matrixPerm->setPermeability( 1., 4 );
    matrixPerm->setPermeability( 1., 8 );
    fracturesPerm->setPermeability( 1.e3, 0 );
    fracturesPerm->setPermeability( 1.e3, 4 );
    fracturesPerm->setPermeability( 1.e3, 8 );
//   const Real kf = 1.e-3; 
//   fracturesPerm->setPermeability( kf, 0 );
    const Real aperture = 1.e-2;
    const Real matrixPoro = 0.25;
    const Real fracturesPoro = 1; 
    propMap.setPropertiesOnMatrix(mesh, matrixPoro, matrixPerm);
    propMap.setPropertiesOnFractures(mesh, aperture, fracturesPoro, fracturesPerm);
    */
    // New Version: data read from file into DataPtr and passed throgh call to the relevant methods.
    propMap.setPropertiesOnMatrix(mesh,dataPtr->getMatrixOtherData(), dataPtr->getMatrixPermeabilityData());
    propMap.setPropertiesOnFractures(mesh,dataPtr->getFractureOtherData(), dataPtr->getFracturePermeabilityData());
//
    std::cout << " done." << std::endl << std::endl;

    std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;

    ExporterVTU exporter;

#ifdef FVCODE3D_EXPORT
    std::cout << "Export..." << std::flush;
    exporter.exportMesh(mesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_mesh.vtu");
    exporter.exportFractures(mesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_fractures.vtu");
    exporter.exportMeshWithFractures(mesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_mesh_fracture.vtu");
    exporter.exportWithProperties(mesh, propMap, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_prop.vtu");
    exporter.exportWireframe(mesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_wire.vtu");
    exporter.exportTetrahedralMesh(mesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_tet.vtu");
    std::cout << " done." << std::endl << std::endl;

    std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;
#endif // FVCODE3D_EXPORT

    std::cout << "Add BCs..." << std::flush;

    auto borders = functions.getBCForStandardDomain();
    /*
    BoundaryConditions::BorderBC leftBC (BorderLabel::Left, Dirichlet, fOne );
    BoundaryConditions::BorderBC rightBC(BorderLabel::Right, Dirichlet, fZero );
    BoundaryConditions::BorderBC backBC (BorderLabel::Back, Neumann, fZero );
    BoundaryConditions::BorderBC frontBC(BorderLabel::Front, Neumann, fZero );
    BoundaryConditions::BorderBC upBC   (BorderLabel::Top, Neumann, fZero );
    BoundaryConditions::BorderBC downBC (BorderLabel::Bottom, Neumann, fZero );

    std::vector<BoundaryConditions::BorderBC> borders;

    borders.push_back( backBC );
    borders.push_back( frontBC );
    borders.push_back( leftBC );
    borders.push_back( rightBC );
    borders.push_back( upBC );
    borders.push_back( downBC );
*/
    BoundaryConditions BC(borders);

    std::cout << " done." << std::endl << std::endl;

    std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


    std::cout << "Assemble rigid mesh..." << std::flush;
    Rigid_Mesh myrmesh(mesh, propMap);
    myrmesh.BuildEdges();
    myrmesh.BuildOrientationFacets();
    std::cout << " done." << std::endl << std::endl;
    

    myrmesh.showMe();
    std::cout << "hMin: " << myrmesh.getMinEdgeSize() << std::endl;
    std::cout << "hMax: " << myrmesh.getMaxEdgeSize() << std::endl;
    std::cout << "hAve: " << myrmesh.getAveEdgeSize() << std::endl;
	std::cout << std::endl;
    std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;

#ifdef FVCODE3D_EXPORT
    std::cout << "Export Fracture Junctures & Tips..." << std::flush;
    exporter.exportFractureJunctures(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_junc.vtu");
    exporter.exportFractureTips(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_tips.vtu");
    exporter.exportEdges(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_edges.vtu");
    exporter.exportFacets(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_facets.vtu");
    std::cout << " done." << std::endl << std::endl;

    std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;
#endif // FVCODE3D_EXPORT

    auto sourceFunction = functions.getSourceFunction();
    std::cout << "Define the problem..." << std::endl<<std::endl;
    
    Pb * darcy(nullptr);
    
		if(dataPtr->getProblemType() == Data::ProblemType::steady && dataPtr->getNumericalMethodType() == Data::NumericalMethodType::MFD)
		{
			darcy = new DarcyPb(dataPtr->getSolverType(), myrmesh, BC, sourceFunction, dataPtr);
		}
		if(dataPtr->getProblemType() == Data::ProblemType::pseudoSteady && dataPtr->getNumericalMethodType() == Data::NumericalMethodType::MFD)
		{
			std::cout<<"Pseudo steady not supported for MFD"<<std::endl<<std::endl;
			return 0;
		}
		if(dataPtr->getNumericalMethodType() == Data::NumericalMethodType::FV && dataPtr->getSolverPolicy() == Data::SolverPolicy::Direct)
		{
			if(dataPtr->getProblemType() == Data::ProblemType::steady)
				darcy = new DarcyPb(dataPtr->getSolverType(), myrmesh, BC, sourceFunction, dataPtr);
			else
				darcy = new PseudoDarcyPb(dataPtr->getSolverType(), myrmesh, BC, sourceFunction, dataPtr);
		}
		if(dataPtr->getNumericalMethodType() == Data::NumericalMethodType::FV && dataPtr->getSolverPolicy() == Data::SolverPolicy::Iterative)
		{
			std::cout<<"Iterative solver not supported for FV"<<std::endl<<std::endl;
			return 0;
		}
			
	
    std::cout << " done." << std::endl << std::endl;

    Real Tsolving1 = 0;

    if(dataPtr->getProblemType() == Data::ProblemType::steady)
    {
		std::cout << "Assemble the problem..." << std::endl<<std::endl;
        darcy->assemble();
        std::cout<<"Assembling done."<<std::endl<<std::endl;
        Tsolving1 = chrono.partial();
        std::cout << "Passed seconds: " << Tsolving1 << " s." << std::endl << std::endl;
        
       if(dataPtr->pressuresInFractures() && dataPtr->getSolverPolicy() == Data::SolverPolicy::Direct)
        {
            FixPressureDofs<DarcyPb> fpd(dynamic_cast<DarcyPb *>(darcy));
            fpd.apply(dataPtr->getPressuresInFractures());
        }
       if(dataPtr->pressuresInFractures() && dataPtr->getSolverPolicy() == Data::SolverPolicy::Iterative)
       {
			std::cout<<"Fix pressure not supported with iterative schemes" <<std::endl;
			return 0;
       }
       
	   if(dynamic_cast<IterativeSolver*>(darcy->getSolverPtr()))
       {
			dynamic_cast<IterativeSolver*>(darcy->getSolverPtr())->setMaxIter(dataPtr->getIterativeSolverMaxIter());
			dynamic_cast<IterativeSolver*>(darcy->getSolverPtr())->setTolerance(dataPtr->getIterativeSolverTolerance());
			dynamic_cast<IterativeSolver*>(darcy->getSolverPtr())->set_precon(dataPtr->getpreconType());
		}

        std::cout << "Solve the problem..." << std::endl<<std::endl;
        darcy->solve();
        
        if(dynamic_cast<IterativeSolver*>(darcy->getSolverPtr()))
        {
              dynamic_cast<IterativeSolver*>(darcy->getSolverPtr())->print();
              std::cout<<std::endl;
        }
    }
    else if(dataPtr->getProblemType() == Data::ProblemType::pseudoSteady)
    {
        dynamic_cast<PseudoDarcyPb *>(darcy)->initialize();
        UInt iter=0;

        std::cout << std::endl << " ... initial solution, t = " << dataPtr->getInitialTime() << " ..." << std::endl;

#ifdef FVCODE3D_EXPORT
        std::stringstream ss;
        ss << dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_solution_";
        ss << iter;
        ss << ".vtu";
        exporter.exportSolution(myrmesh, ss.str(), dynamic_cast<PseudoDarcyPb *>(darcy)->getOldSolution());
        ss.str(std::string());
        ss << dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_solution_f_";
        ss << iter;
        ss << ".vtu";
        exporter.exportSolutionOnFractures(myrmesh, ss.str(), dynamic_cast<PseudoDarcyPb *>(darcy)->getOldSolution());
        std::cout << "  done." << std::endl << std::endl;
#endif // FVCODE3D_EXPORT

        ++iter;

        for(Real t = dataPtr->getInitialTime() + dataPtr->getTimeStep() ; t <= dataPtr->getEndTime(); t+=dataPtr->getTimeStep(), ++iter)
        {
            std::cout << " ... at t = " << t << " ..." << std::endl;
            darcy->assemble();
            if(dataPtr->pressuresInFractures())
            {
                FixPressureDofs<PseudoDarcyPb> fpd(dynamic_cast<PseudoDarcyPb *>(darcy));
                fpd.apply(dataPtr->getPressuresInFractures());
            }

            darcy->solve();

            if(dynamic_cast<IterativeSolver*>(&(darcy->getSolver())))
            {
                std::cout << std::endl;
                std::cout << "\t# iterations: " << dynamic_cast<IterativeSolver*>(darcy->getSolverPtr())->getIter() << std::endl;
                std::cout << "\tResidual: " << dynamic_cast<IterativeSolver*>(darcy->getSolverPtr())->getResidual() << std::endl;
                std::cout << std::endl;
            }

#ifdef FVCODE3D_EXPORT
            ss.str(std::string());
            std::cout << " Export Solution" << std::flush;
            ss << dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_solution_";
            ss << iter;
            ss << ".vtu";
            exporter.exportSolution(myrmesh, ss.str(), darcy->getSolver().getSolution());
            ss.str(std::string());
            ss << dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_solution_f_";
            ss << iter;
            ss << ".vtu";
            exporter.exportSolutionOnFractures(myrmesh, ss.str(), darcy->getSolver().getSolution());
            std::cout << "  done." << std::endl << std::endl;
#endif // FVCODE3D_EXPORT
        }
    }
	
	std::cout<< "Time to solve the system: "<< chrono.partial()-Tsolving1 <<"s."<< std::endl << std::endl;
    std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;  
     
    if(dataPtr->getNumericalMethodType() == Data::MFD){
		//UInt numFacetsTot   = myrmesh.getFacetsVector().size() + myrmesh.getFractureFacetsIdsVector().size();
		UInt numCellsTot    = myrmesh.getCellsVector().size() + myrmesh.getFractureFacetsIdsVector().size();
		std::cout << "Export Solution..." << std::flush;
		exporter.exportSolution( myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_solution.vtu", 
			darcy->getSolver().getSolution().tail(numCellsTot) );
		std::cout << "Export Solution on Fractures..." << std::flush;
		exporter.exportSolutionOnFractures( myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_solution_f.vtu", 
			darcy->getSolver().getSolution().tail(numCellsTot) );
		std::cout << " done." << std::endl << std::endl;
	}
	    
    else{
		std::cout << "Export Solution..." << std::flush;
		exporter.exportSolution(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_solution.vtu", darcy->getSolver().getSolution());
		std::cout << "Export Solution on Fractures..." << std::flush;
		exporter.exportSolutionOnFractures(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_solution_f.vtu",
        darcy->getSolver().getSolution());
        std::cout << " done." << std::endl << std::endl;
	}

    std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;

    if(dataPtr->getDumpMatrices())
      {
        std::cout<<"  *** DUMPING SYSTEM MATRICES ***"<<std::endl;
        auto solverSPtr=darcy->getSolverSharedPtr();
        auto solverPtr=dynamic_cast<IterativeSolver*>(solverSPtr.get());
        if(solverPtr == nullptr)
          {
            std::cerr<<" CANNOT DUMP MATRIX FOR DIRECT SOLVER, NOT YET IMPLEMENTED ***"<<std::endl;
          }
        else
          {
            bool status=dumpSaddlePointMatrix(solverPtr->getA());
            if (!status)
              {
                std::cerr<<"SOMETHING WRONG DUMPING SYSTEM MATRICES"<<std::endl;
              }
            else
              {
                std::cout<<"Dumped M, B and T"<<std::endl;
              }
            bool lumped = dataPtr->getLumpedMimetic();
            status=dumpApproximateSchurMatrix(solverPtr->getA(),lumped);
            if (!status)
              {
                std::cerr<<"SOMETHING WRONG DUMPING SYSTEM MATRICES"<<std::endl;
              }
            else
              {
                std::cout<<"Dumped Approximate M and approximate Schur"<<std::endl;
              }

          }
      }


#ifdef FVCODE3D_EXPORT
    std::cout << "Export Property..." << std::flush;
    exporter.exportWithProperties(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_APP.vtu", Aperture | Permeability | Porosity);
    std::cout << " done." << std::endl << std::endl;

    std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;
#endif // FVCODE3D_EXPORT

    delete darcy;
    delete importer;

    chrono.stop();
    return 0;
}

