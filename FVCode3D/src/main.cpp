/*! @mainpage
 *
 * @section intro Introduction
 *
 * This code allows to solve the single-phase flow in a fractured porous medium modeled via the Darcy equation by using the
 Finite Volume (FV) method with a two-point flux approximation (TPFA). The fractures are required to be matching with the grid.\n
 * The code is written in C++.
 * As input, the program requires a mesh as '.fvg' format file.
 * The output is given as '.vtu' file format, viewable by Paraview.\n
 *
 * @section howto How to execute the program
 *
 * To run the program just execute\n
 * @code
 * ./fvcode3d.exe
 * @endcode
 * By default, the parameters need are read from data.txt file.\n
 * With the option '-f ', it is possible to specify the data file.\n
 * Ex.\n
 * @code
 * ./fvcode3d.exe -f data.txt
 * @endcode
 *
 * @section input Input parameters
 *
 * The program requires:\n
 * - an input mesh in '.fvg' format\n
 * - the data file\n
 *
 * The file format '.fvg' consists of four parts:\n
 * - list of POINTS, each of them described by its coordinates x, y, z;\n
 * - list of FACETS, each of them described by a list of ids of POINTS, plus some properties if the facet represents a
 fracture;\n
 * - list of CELLS, each of them described by a list of ids of FACETS, plus some properties;\n
 * - list of FRACTURE NETWORKS, each of them described by a list of FACETS that belong to the same network.\n
 *
 * The properties related to the facets that represent a fracture are the permeability, the porosity and the aperture.\n
 * The properties related to the cells (porous medium) are the permeability and the porosity.
 *
 * The data file consists of six sections:
 * - @c mesh: the input file parameters such as mesh location and filename\n
 * - @c output: the output directory and filename\n
 * - @c problem: parameters related to the problem, such as the problem type (steady or unsteady) or if apply or not the
 source/sink term\n
 * - @c fluid: mobility and compressibility of the fluid\n
 * - @c bc: the rotation to apply around the z axis, needed to correctly apply the boundary conditions\n
 * - @c solver: the solver to use and some parameters such as the tolerance and the number of iterations\n
 *
 * @section class Main Classes
 *
 * The main classes involved are the following:\n
 *
 * - @a Data: collects all the parameters used by the program that can be set through the datafile;\n
 * - @a Importer: reads the input file mesh;\n
 * - @a PropertiesMap: collects the properties of the fractures and porous medium;\n
 * - @a Mesh3D: stores the geometrical mesh as @a Point3D, @a Facet3D and @a Cell3D;\n
 * - @a Rigid_Mesh: converts a @a Mesh3D in a rigid format, suitable for assemble the problem;\n
 * - @a BoundaryConditions: defines the boundary conditions to assign to the problem; BCs can be Dirichlet or Neumann;\n
 * - @a Problem: defines the problem to solve: it can be a steady problem (@a DarcySteady) or a time-dependent problem
 (@a DarcyPseudoSteady); it assembles the matrix and the right hand side;\n
 * - @a Solver: solves the @a Problem; the code implements direct solvers as well as iterative solvers;\n
 * - @a Exporter: exports as '.vtu' files the mesh, the properties, the solution and more.\n
 *
 * @section es Example of usage
 *
 * Here, we briefly describe how to use the code.\n
 * First of all, read the necessary parameters from the datafile:\n
 *
 * @code
 *
 * // we use the GetPot utility to access the datafile
 * GetPot command_line(argc,argv);
 * const std::string dataFileName = command_line.follow("data.txt", 2, "-f", "--file");
 * // Read Data from the datafile
 * DataPtr_Type dataPtr(new Data(dataFileName));
 *
 * @endcode
 *
 * Then we define the Mesh3D, the PropertiesMap so that we can import the mesh:
 *
 * @code
 *
 * //Define Mesh and Properties
 * Mesh3D mesh;
 * PropertiesMap propMap(dataPtr->getMobility(), dataPtr->getCompressibility());
 * //Create Importer
 * Importer * importer = 0;
 * importer = new ImporterForSolver(dataPtr->getMeshDir() + dataPtr->getMeshFile(), mesh, propMap);
 * //Import grid file
 * importer->import(dataPtr->fractureOn());
 *
 * @endcode
 *
 * After the reading of the mesh, it is necessary to perform some operations to process the mesh:
 *
 * @code
 *
 * //Compute the cells that separate each facet
 * mesh.updateFacetsWithCells();
 * //Compute the neighboring cells of each cell
 * mesh.updateCellsWithNeighbors();
 * //Set labels on boundary (necessary for BCs)
 * importer->extractBC(dataPtr->getTheta());
 * //Compute facet ids of the fractures (creates fracture networks)
 * mesh.updateFacetsWithFractures();
 *
 * @endcode
 *
 * We can now create the boundary conditions:
 *
 * @code
 *
 * //Assign to the border marked as 'Left' a Dirichlet condition equal to 1
 * BoundaryConditions::BorderBC leftBC (BorderLabel::Left, Dirichlet, fOne );
 * //And so on...
 * BoundaryConditions::BorderBC rightBC(BorderLabel::Right, Dirichlet, fZero );
 * ...
 *
 * //We create a vector that contains all the BorderBC...
 * std::vector<BoundaryConditions::BorderBC> borders;
 * //...and we insert into it the BCs previously created
 * borders.push_back( backBC );
 * borders.push_back( frontBC );
 * ...
 *
 * //Finally we create the BoundaryConditions
 * BoundaryConditions BC(borders);
 *
 * @endcode
 *
 * Then we create a Rigid_Mesh
 * @code
 *
 * //The Rigid_Mesh requires the Mesh3D and the PropertesMap
 * Rigid_Mesh myrmesh(mesh, propMap);
 *
 * @endcode
 *
 * We can proceed by building the problem:
 *
 * @code
 *
 * Problem<CentroidQuadrature, CentroidQuadrature> * darcy(nullptr);
 * typedef DarcySteady<CentroidQuadrature, CentroidQuadrature> DarcyPb;
 * typedef DarcyPseudoSteady<CentroidQuadrature, CentroidQuadrature, SpMat, TimeScheme::BDF2> PseudoDarcyPb;
 *
 * //If the problem is steady
 * if(dataPtr->getProblemType() == Data::ProblemType::steady)
 * {
 *     darcy = new DarcyPb(dataPtr->getSolverType(), myrmesh, BC, SS, dataPtr);
 * }
 * //If the problem is unsteady
 * else if(dataPtr->getProblemType() == Data::ProblemType::pseudoSteady)
 * {
 *     darcy = new PseudoDarcyPb(dataPtr->getSolverType(), myrmesh, BC, SS, dataPtr);
 * }
 *
 * //If the solver used is an iterative one, set tolerance and the maximum number of iterations
 * if(dynamic_cast<IterativeSolver*>(darcy->getSolverPtr()))
 * {
 *     dynamic_cast<IterativeSolver*>(darcy->getSolverPtr())->setMaxIter(dataPtr->getIterativeSolverMaxIter());
 *     dynamic_cast<IterativeSolver*>(darcy->getSolverPtr())->setTolerance(dataPtr->getIterativeSolverTolerance());
 * }
 *
 * @endcode
 *
 * Now, we can solve the problem:
 *
 * @code
 *
 * //In the case of a steady problem, we assemble and the solve the problem
 * if(dataPtr->getProblemType() == Data::ProblemType::steady)
 * {
 *     darcy->assemble();
 *     if(dataPtr->pressuresInFractures())
 *     {
 *         FixPressureDofs<DarcyPb> fpd(dynamic_cast<DarcyPb *>(darcy));
 *         fpd.apply(dataPtr->getPressuresInFractures());
 *     }
 *     darcy->solve();
 * }
 * //In the case of a unsteady problem, we initialize the problem and for each time step, we assemble and solve the
 problem
 * else if(dataPtr->getProblemType() == Data::ProblemType::pseudoSteady)
 * {
 *     dynamic_cast<PseudoDarcyPb *>(darcy)->initialize();
 *     for(Real t = dataPtr->getInitialTime() + dataPtr->getTimeStep() ; t <= dataPtr->getEndTime(); t+=dataPtr->getTimeStep())
 *     {
 *         darcy->assemble();
 *         if(dataPtr->pressuresInFractures())
 *         {
 *             FixPressureDofs<PseudoDarcyPb> fpd(dynamic_cast<PseudoDarcyPb *>(darcy));
 *             fpd.apply(dataPtr->getPressuresInFractures());
 *         }
 *         darcy->solve();
 *     }
 * }
 *
 * @endcode
 *
 * The solution can be exported:
 *
 * @code
 *
 * //Define Exporter
 * ExporterVTU exporter;
 * //Export solution on matrix and on the fractures
 * exporter.exportSolution(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_solution.vtu", darcy->getSolver().getSolution());
 * exporter.exportSolutionOnFractures(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_solution_f.vtu", darcy->getSolver().getSolution());
 *
 * @endcode
 *
 * Remember to delete the Problem and the Importer instance:
 *
 * @code
 *
 * delete darcy;
 * delete importer;
 *
 * @endcode
 *
 */

//#define NDEBUG // add this macro to disable asserts
#include <cassert>
#include <iostream>
#include <chrono>
#include <vector>
#include <fstream>

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

//    std::cout << "Set uniform properties..." << std::flush;
//    std::shared_ptr<PermeabilityBase> matrixPerm( new PermeabilityDiagonal );
//    std::shared_ptr<PermeabilityBase> fracturesPerm( new PermeabilityScalar );
//    matrixPerm->setPermeability( 1., 0 );
//    matrixPerm->setPermeability( 1., 4 );
//    matrixPerm->setPermeability( 1., 8 );
//    fracturesPerm->setPermeability( 1.e6, 0 );
//    propMap.setPropertiesOnMatrix(mesh, 0.25, matrixPerm);
//    propMap.setPropertiesOnFractures(mesh, 1e-2, 1, fracturesPerm);
//    std::cout << " done." << std::endl << std::endl;

    std::cout << "Export..." << std::flush;
    ExporterVTU exporter;
    exporter.exportMesh(mesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_mesh.vtu");
    exporter.exportFractures(mesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_fractures.vtu");
    exporter.exportMeshWithFractures(mesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_mesh_fracture.vtu");
    exporter.exportWithProperties(mesh, propMap, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_prop.vtu");
    exporter.exportWireframe(mesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_wire.vtu");
    exporter.exportTetrahedralMesh(mesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_tet.vtu");
//  if(dataPtr->getMeshType() == Data::MeshFormatType::TPFA)
//      saveAsSolverFormat(dataPtr->getMeshDir() + dataPtr->getMeshFile() + "_new.fvg", mesh, propMap);
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

    std::cout << "Solve problem..." << std::flush;
    if(dataPtr->getProblemType() == Data::ProblemType::steady)
    {
        darcy->assemble();
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
        darcy->solve();
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
                std::cout << "\t# iterations: " << dynamic_cast<IterativeSolver*>(darcy->getSolverPtr())->getIter() << std::endl;
                std::cout << "\tResidual: " << dynamic_cast<IterativeSolver*>(darcy->getSolverPtr())->getResidual() << std::endl;
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


    std::cout << "Export Solution..." << std::flush;
    exporter.exportSolution(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_solution.vtu", darcy->getSolver().getSolution());
    std::cout << " done." << std::endl << std::endl;

    std::cout << "Export Solution on Fractures..." << std::flush;
    exporter.exportSolutionOnFractures(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_solution_f.vtu", darcy->getSolver().getSolution());
    std::cout << " done." << std::endl << std::endl;

    if(dataPtr->MSROn())
    {
        std::cout << "Export SubRegions..." << std::flush;
        exporter.exportSolution(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_colour.vtu",
    multipleSubRegions->getColorsVector() );
        std::cout << " done." << std::endl << std::endl;
    }

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
