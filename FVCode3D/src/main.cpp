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
#include <FVCode3D/multipleSubRegions/MultipleSubRegions.hpp>
#include <FVCode3D/utility/Evaluate.hpp>
#include "functions.hpp"

using namespace FVCode3D;

typedef Problem<CentroidQuadrature, CentroidQuadrature> Pb;
typedef DarcySteady<CentroidQuadrature, CentroidQuadrature> DarcyPb;
typedef DarcyPseudoSteady<CentroidQuadrature, CentroidQuadrature, SpMat, TimeScheme::BDF2> PseudoDarcyPb;

int main(int argc, char * argv[])
{
    GetPot command_line(argc,argv);
    const std::string dataFileName = command_line.follow("data.txt", 2, "-f", "--file");

    Chrono chrono;
    chrono.start();

    std::cout << "Read Data from " << dataFileName << "..." << std::flush;
    DataPtr_Type dataPtr(new Data(dataFileName));
    std::cout << " done." << std::endl;

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

    std::cout << "Set uniform properties..." << std::flush;
    std::shared_ptr<PermeabilityBase> matrixPerm( new PermeabilityDiagonal );
    std::shared_ptr<PermeabilityBase> fracturesPerm( new PermeabilityScalar );
    matrixPerm->setPermeability( 1., 0 );
    matrixPerm->setPermeability( 1., 4 );
    matrixPerm->setPermeability( 1., 8 );
    const Real kf = 1.e3 ; 
    fracturesPerm->setPermeability( kf, 0 );
    const Real aperture = 1.e-2;
    const Real matrixPoro = 0.25;
    const Real fracturesPoro = 1; 
    propMap.setPropertiesOnMatrix(mesh, matrixPoro, matrixPerm);
    propMap.setPropertiesOnFractures(mesh, aperture, fracturesPoro, fracturesPerm);
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

    BoundaryConditions BC(borders);
    std::cout << " done." << std::endl << std::endl;

    std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


    std::cout << "Assemble rigid mesh..." << std::flush;
    Rigid_Mesh myrmesh(mesh, propMap);
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

    std::cout << "Define the problem..." << std::flush;
    Pb * darcy(nullptr);
    if(dataPtr->getProblemType() == Data::ProblemType::steady)
    {
        darcy = new DarcyPb(dataPtr->getSolverType(), myrmesh, BC, SS, dataPtr);
    }
    else if(dataPtr->getProblemType() == Data::ProblemType::pseudoSteady)
    {
        darcy = new PseudoDarcyPb(dataPtr->getSolverType(), myrmesh, BC, SS, dataPtr);
    }
    std::cout << " done." << std::endl << std::endl;

    if(dynamic_cast<IterativeSolver*>(darcy->getSolverPtr()))
    {
        dynamic_cast<IterativeSolver*>(darcy->getSolverPtr())->setMaxIter(dataPtr->getIterativeSolverMaxIter());
        dynamic_cast<IterativeSolver*>(darcy->getSolverPtr())->setTolerance(dataPtr->getIterativeSolverTolerance());
    }

    MSR<Pb> * multipleSubRegions(nullptr);

    if(dataPtr->getProblemType() == Data::ProblemType::steady)
    {
		std::cout << "Assemble the problem..." << std::endl<<std::endl;
        darcy->assemble();
        std::cout<<"Assembling done."<<std::endl<<std::endl;
        std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;
        
        if(dataPtr->MSROn())
        {
            multipleSubRegions = new MSR<Pb>(darcy, dataPtr);
            multipleSubRegions->setup();
        }
        if(dataPtr->pressuresInFractures())
        {
            FixPressureDofs<DarcyPb> fpd(dynamic_cast<DarcyPb *>(darcy));
            fpd.apply(dataPtr->getPressuresInFractures());
        }
        std::cout << "Solve the problem..." << std::endl<<std::endl;
        darcy->solve();
        std::cout<<"done."<<std::endl<<std::endl;
        
        if(dynamic_cast<IterativeSolver*>(darcy->getSolverPtr()))
        {
            std::cout << std::endl;
            std::cout << "\t# iterations: " << dynamic_cast<IterativeSolver*>(darcy->getSolverPtr())->getIter() << std::endl;
            std::cout << "\tResidual: " << dynamic_cast<IterativeSolver*>(darcy->getSolverPtr())->getResidual() << std::endl;
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

    if(dataPtr->MSROn())
    {
        std::cout << "Compute sub-regions..." << std::flush;
        multipleSubRegions->createRegions();
        std::cout << " done." << std::endl;
        std::cout << "Compute transmissibilities..." << std::flush;
        multipleSubRegions->computeTransmissibility();
        std::cout << " done." << std::endl << std::endl;
    }

    std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;
    
    
    if(dataPtr->getNumericalMethodType() == Data::MFD){
		UInt numFacetsTot   = myrmesh.getFacetsVector().size() + myrmesh.getFractureFacetsIdsVector().size();
		UInt numCellsTot    = myrmesh.getCellsVector().size() + myrmesh.getFractureFacetsIdsVector().size();
		std::cout << "Export Solution..." << std::flush;
		exporter.exportSolution( myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_solution.vtu", 
			darcy->getSolver().getSolution().segment(numFacetsTot,numCellsTot) );
		std::cout << "Export Solution on Fractures..." << std::flush;
		exporter.exportSolutionOnFractures( myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_solution_f.vtu", 
			darcy->getSolver().getSolution().segment(numFacetsTot,numCellsTot) );
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

#ifdef FVCODE3D_EXPORT
    if(dataPtr->MSROn())
    {
        std::cout << "Export SubRegions..." << std::flush;
        exporter.exportSolution(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_colour.vtu",
    multipleSubRegions->getColorsVector() );
        std::cout << " done." << std::endl << std::endl;
    }
#endif // FVCODE3D_EXPORT

    std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;

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
