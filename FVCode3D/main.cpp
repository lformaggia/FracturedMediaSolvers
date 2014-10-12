//#define NDEBUG // add this macro to disable asserts
//#define FluxFracture
#include <cassert>
#include <iostream>
#include <chrono>
#include <vector>
#include <fstream>

#include <unsupported/Eigen/SparseExtra>

#include "core/TypeDefinition.hpp"
#include "core/data.hpp"
//#include "geometry/"
#include "mesh/Rigid_Mesh.hpp"
#include "property/Properties.hpp"
#include "mesh/cartesianGrid.hpp"
#include "boundaryCondition/BC.hpp"
#include "quadrature/Quadrature.hpp"
#include "import/import.hpp"
#include "export/exportVTU.hpp"
#include "utility/converter.hpp"
#include "solver/solver.hpp"
#include "problem/problem.hpp"
#include "problem/darcySteady.hpp"
#include "problem/darcyPseudoSteady.hpp"
#include "assembler/fixPressureDofs.hpp"
#include "functions.hpp"

typedef Problem<EigenBiCGSTAB, CentroidQuadrature, CentroidQuadrature> Pb;
typedef DarcySteady<EigenBiCGSTAB, CentroidQuadrature, CentroidQuadrature> DarcyPb;
typedef DarcyPseudoSteady<EigenBiCGSTAB, CentroidQuadrature, CentroidQuadrature, TimeScheme::BDF2> PseudoDarcyPb;

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
    Geometry::Mesh3D mesh;
    Geometry::PropertiesMap propMap(dataPtr->getMobility(), dataPtr->getCompressibility());
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
    BoundaryConditions::BorderBC backBC (1, Neumann, fMinusOne );
    BoundaryConditions::BorderBC frontBC(2, Dirichlet, fOne );
    BoundaryConditions::BorderBC leftBC (3, Neumann, fZero );
    BoundaryConditions::BorderBC rightBC(4, Neumann, fZero );
    BoundaryConditions::BorderBC upBC   (5, Neumann, fZero );
    BoundaryConditions::BorderBC downBC (6, Neumann, fZero );

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


    std::cout << "Export Fracture Junctures & Tips..." << std::flush;
    exporter.exportFractureJunctures(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_junc.vtu");
    exporter.exportFractureTips(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_tips.vtu");
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


    // Save matrix and rhs
    Eigen::saveMarket(darcy->getA(),"stiffness.m");
    Eigen::saveMarket(darcy->getb(),"rhs.m");

    // Transmissibility fracture-frature and boundary fluxes
    const std::vector<Geometry::Point3D> & meshNodes = myrmesh.getNodesVector();
    const Vector & solution = darcy->getSolver().getSolution();
    const UInt nbCells = myrmesh.getCellsVector().size();
    Real pressure1(0.), pressure2(0.);
    Real fracture1Area(0.), fracture2Area(0.);
    Real fractureFlux(0.), inletFlux(0.), outletFlux(0.), fractureFluxIn(0.), fractureFluxOut(0.);
    Real h( std::numeric_limits<Real>::max() );
    const Real midX(1.);

    Darcy::StiffMatrix StiffM(myrmesh, BC);
    StiffM.assemble();
    StiffM.reconstructFlux(solution);
    const Vector & flux = StiffM.getFlux();

    std::cout << " Export Flux" << std::flush;
    exporter.exportFlux(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_flux.vtu", flux);
    //exporter.exportFluxOnFractures(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_flux_f.vtu", flux);
    std::cout << "  done." << std::endl << std::endl;

    // Boundary fluxes
    for ( auto facet_it : myrmesh.getBorderFacetsIdsVector() )
    {
        const UInt borderID = facet_it.getBorderId();
        const UInt facetID = facet_it.getFacetId();
        if (borderID == 1 || borderID == 2)
        {
            const UInt cellId = facet_it.getSeparated()[0];
            const Real pD = BC.getBordersBCMap().at(borderID).getBC()(facet_it.getCentroid());
            const Real pC = solution[cellId];
            const Geometry::Point3D UNormal = facet_it.getUNormal();
            const Geometry::Point3D mGradP = (pD > pC) ?
                                (myrmesh.getCellsVector()[cellId].getCentroid() - facet_it.getCentroid()) :
                                (facet_it.getCentroid() - myrmesh.getCellsVector()[cellId].getCentroid());
            const Real alpha = (dotProduct(UNormal, mGradP) >= 0.) ? 1. : -1. ;
            const Real fluxC = alpha * flux[facetID];

            if(borderID == 1)
            {
                inletFlux += fluxC;
            }
            if(borderID == 2)
            {
                outletFlux += fluxC;
            }

        }
    }

    const Real deltaFlux( std::abs(inletFlux) - std::abs(outletFlux) );

    std::cout << std::setprecision(16)
              << " inletFlux " << std::abs(inletFlux)
              << std::endl;

    std::cout << std::setprecision(16)
              << " outletFlux " << std::abs(outletFlux)
              << std::endl;

    std::cout << std::setprecision(16)
              << " deltaFlux " << std::abs(deltaFlux)
              << std::endl;

    // Transmissibility fracture-frature
#ifdef FluxFracture
    for ( auto facet_it : myrmesh.getFractureFacetsIdsVector() )
    {
        const auto& facetNeighbors = facet_it.getFractureNeighbors ();
        const UInt facetIdasCell = facet_it.getIdasCell();
        const Geometry::Point3D facetCentre = facet_it.getCentroid();
        const Real facetArea = facet_it.getSize();

        h = std::min( h, std::sqrt( facetArea ) / 2. );

        if ( facetCentre[0] < midX )
        {
            pressure1 += solution[ facetIdasCell ] * facetArea;
            fracture1Area += facetArea;
        } // if

        if ( facetCentre[0] > midX /*&& facetCentre[0] < 4100.0*/ )
        {
            pressure2 += solution[ facetIdasCell ] * facetArea;
            fracture2Area += facetArea;
        } // if

            const Real F_permeability = propMap.getProperties(facet_it.getZoneCode()).M_permeability;
            const Real F_aperture = propMap.getProperties(facet_it.getZoneCode()).M_aperture;
            const Real Df = F_aperture/2.;
            const Real alphaf = facet_it.getSize() * F_permeability / Df;

            const Real alpha1 = StiffM.Findalpha(facet_it.getSeparated()[0], &facet_it);
            const Real alpha2 = StiffM.Findalpha(facet_it.getSeparated()[1], &facet_it);

            const Real T1f = alpha1*alphaf/(alpha1 + alphaf) * propMap.getMobility();
            const Real flux1 = T1f * (solution[facet_it.getSeparated()[0]] - solution[facetIdasCell]) ;

            const Real T2f = alphaf*alpha2/(alphaf + alpha2) * propMap.getMobility();
            const Real flux2 = T2f * (solution[facet_it.getSeparated()[1]] - solution[facetIdasCell]) ;

            if (flux1 > 0)
                fractureFluxIn += flux1;
            else
                fractureFluxOut += flux1;

            if (flux2 > 0)
                fractureFluxIn += flux2;
            else
                fractureFluxOut += flux2;

        for ( auto& facetNeighbors_it : facetNeighbors )
        {
            const UInt node1Id = facetNeighbors_it.first.first;
            const UInt node2Id = facetNeighbors_it.first.second;

            if ( meshNodes[node1Id][0] == midX && meshNodes[node2Id][0] == midX && facetCentre[0] < midX )
            {
                const UInt neighborFacetIdasCell = facetNeighbors_it.second[0] + nbCells;
                const Real trasmissibility = darcy->getA().coeffRef( facetIdasCell, neighborFacetIdasCell );
                const Real deltaSolution = solution[ facetIdasCell ] - solution[ neighborFacetIdasCell ];
                fractureFlux += trasmissibility * deltaSolution;
            } // if
        } // for
    } // for

    std::cout << "Fracture mesh size " << h << std::endl;

    const Real deltaFractureFlux( std::abs(fractureFluxIn) - std::abs(fractureFluxOut) );

    std::cout << std::setprecision(16)
                << " inletFractureFlux " << std::abs(fractureFluxIn)
                << std::endl;

    std::cout << std::setprecision(16)
                << " outletFractureFlux " << std::abs(fractureFluxOut)
                << std::endl;

    std::cout << std::setprecision(16)
                << " deltaFractureFlux " << std::abs(deltaFractureFlux)
                << std::endl;

    pressure1 /= fracture1Area;
    pressure2 /= fracture2Area;

    std::cout << std::setprecision(16)
                << "size1 " << fracture1Area
                << std::endl
                << " size2 " << fracture2Area
                << std::endl;

    const Real deltaFracturePressure( pressure1 - pressure2 );

    std::cout << std::setprecision(16)
                << "pressure1 " << pressure1
                << std::endl
                << " pressure2 " << pressure2
                << " delta " << deltaFracturePressure
                << std::endl;

    const Real fractureTrasmissibility( fractureFlux / deltaFracturePressure );

    const Real T_teo = 1.1250e+04;

    std::cout << std::setprecision(16)
                << " flux " << std::abs(fractureFlux)
                << std::endl;

    std::cout << std::setprecision(16)
                << " trasmissibility " << std::abs(fractureTrasmissibility)
                << std::endl;

    std::cout << std::setprecision(16)
                << " err " << std::abs((T_teo - std::abs(fractureTrasmissibility)) / T_teo)
                << std::endl;

    std::cout << std::setprecision(16)
                << " Et " << 2*1.e-16 * std::max(pressure1, pressure2) / ( (pressure1 - pressure2)*(pressure1 - pressure2))
                << std::endl;
#endif

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
