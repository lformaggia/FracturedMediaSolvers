//#define NDEBUG // add this macro to disable asserts
#include <cassert>
#include <iostream>
#include <chrono>
#include <vector>
#include <fstream>

#define FVCODE3D_HAS_UMFPACK

#include "core/TypeDefinition.hpp"
#include "core/Data.hpp"
#include "mesh/RigidMesh.hpp"
#include "property/Properties.hpp"
#include "mesh/CartesianGrid.hpp"
#include "boundaryCondition/BC.hpp"
#include "quadrature/Quadrature.hpp"
#include "import/Import.hpp"
#include "export/ExportVTU.hpp"
#include "utility/Converter.hpp"
#include "solver/Solver.hpp"
#include "problem/Problem.hpp"
#include "problem/DarcySteady.hpp"
#include "problem/DarcyPseudoSteady.hpp"
#include "assembler/FixPressureDofs.hpp"
//#include "multipleSubRegions/MultipleSubRegions.hpp"
#include "functions.hpp"

using namespace FVCode3D;

typedef EigenLU SolverType;

typedef Problem<SolverType, CentroidQuadrature, CentroidQuadrature> Pb;
typedef DarcySteady<SolverType, CentroidQuadrature, CentroidQuadrature> DarcyPb;
typedef DarcyPseudoSteady<SolverType, CentroidQuadrature, CentroidQuadrature, SpMat, TimeScheme::BDF2> PseudoDarcyPb;

int main(int argc, char * argv[])
{
    GetPot command_line(argc,argv);
    const std::string dataFileName = command_line.follow("data.txt", 2, "-f", "--file");

    Chrono chrono;
    chrono.start();

    std::cout << "Read Data..." << std::flush;
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

    propMap.setPropertiesOnMatrix(mesh, 0.25, 1);
    propMap.setPropertiesOnFractures(mesh, 1e-2, 1, 1e6);

    std::cout << "Export..." << std::flush;
    ExporterVTU exporter;
    exporter.exportMesh(mesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_mesh.vtu");
    exporter.exportFractures(mesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_fractures.vtu");
    exporter.exportMeshWithFractures(mesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_mesh_fracture.vtu");
    exporter.exportWithProperties(mesh, propMap, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_prop.vtu");
    exporter.exportWireframe(mesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_wire.vtu");
    exporter.exportTetrahedralMesh(mesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_tet.vtu");
//  if(dataPtr->getMeshType() == Data::MeshFormatType::TPFA)
//      saveAsSolverFormat(dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_new.fvg", mesh, propMap);
    std::cout << " done." << std::endl << std::endl;

    std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


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

    std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


    std::cout << "Export Fracture Junctures & Tips..." << std::flush;
    exporter.exportFractureJunctures(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_junc.vtu");
    exporter.exportFractureTips(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_tips.vtu");
    exporter.exportEdges(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_edges.vtu");
    exporter.exportFacets(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_facets.vtu");
    std::cout << " done." << std::endl << std::endl;

    std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


    std::cout << "Build problem..." << std::flush;
    Pb * darcy(nullptr);
    if(dataPtr->getProblemType() == Data::ProblemType::steady)
        darcy = new DarcyPb(myrmesh, BC, SS, dataPtr);
    else if(dataPtr->getProblemType() == Data::ProblemType::pseudoSteady)
        darcy = new PseudoDarcyPb(myrmesh, BC, SS, dataPtr);
    std::cout << " done." << std::endl << std::endl;

    if(dynamic_cast<IterativeSolver*>(&(darcy->getSolver())))
    {
        dynamic_cast<IterativeSolver*>(&(darcy->getSolver()))->setMaxIter(1000);
        dynamic_cast<IterativeSolver*>(&(darcy->getSolver()))->setTolerance(1e-8);
    }

    std::cout << "Solve problem..." << std::flush;
    if(dataPtr->getProblemType() == Data::ProblemType::steady)
    {
        darcy->assemble();
        if(dataPtr->pressuresInFractures())
        {
            FixPressureDofs<DarcyPb> fpd(dynamic_cast<DarcyPb *>(darcy));
            fpd.apply(dataPtr->getPressuresInFractures());
        }
        darcy->solve();
        if(dynamic_cast<IterativeSolver*>(&(darcy->getSolver())))
        {
            std::cout << std::endl;
            std::cout << "\t# iterations: " << dynamic_cast<IterativeSolver*>(&(darcy->getSolver()))->getIter() << std::endl;
            std::cout << "\tResidual: " << dynamic_cast<IterativeSolver*>(&(darcy->getSolver()))->getResidual() << std::endl;
        }
    }
    else if(dataPtr->getProblemType() == Data::ProblemType::pseudoSteady)
    {
        dynamic_cast<PseudoDarcyPb *>(darcy)->initialize();
        UInt iter=0;

        std::cout << std::endl << " ... initial solution, t = " << dataPtr->getInitialTime() << " ..." << std::endl;

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
                std::cout << "\t# iterations: " << dynamic_cast<IterativeSolver*>(&(darcy->getSolver()))->getIter() << std::endl;
                std::cout << "\tResidual: " << dynamic_cast<IterativeSolver*>(&(darcy->getSolver()))->getResidual() << std::endl;
                std::cout << std::endl;
            }

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
        }
    }
    std::cout << " done." << std::endl << std::endl;

    std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


    std::cout << "Export Solution..." << std::flush;
    exporter.exportSolution(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_solution.vtu", darcy->getSolver().getSolution());
    std::cout << " done." << std::endl << std::endl;

    std::cout << "Export Solution on Fractures..." << std::flush;
    exporter.exportSolutionOnFractures(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_solution_f.vtu", darcy->getSolver().getSolution());
    std::cout << " done." << std::endl << std::endl;

    std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;


    std::cout << "Export Property..." << std::flush;
    exporter.exportWithProperties(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_APP.vtu", Aperture | Permeability | Porosity);
    std::cout << " done." << std::endl << std::endl;

    std::cout << "Passed seconds: " << chrono.partial() << " s." << std::endl << std::endl;

    delete darcy;
    delete importer;

    chrono.stop();

    return 0;
}
