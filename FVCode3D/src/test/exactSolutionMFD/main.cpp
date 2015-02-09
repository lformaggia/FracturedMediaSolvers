#include <FVCode3D/FVCode3D.hpp>

// Add this macro the enable the exporter
#define FVCODE3D_EXPORT

using namespace FVCode3D;

typedef Problem<CentroidQuadrature, CentroidQuadrature> Pb;
typedef DarcySteady<CentroidQuadrature, CentroidQuadrature> DarcyPb;

int main (int argc, char** argv)
{
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

    std::shared_ptr<PermeabilityBase> matrixPerm( new PermeabilityScalar );
    matrixPerm->setPermeability( 1., 0 );
    const Real matrixPoro = 0.25;
    propMap.setPropertiesOnMatrix(mesh, matrixPoro, matrixPerm);

    Func fZero  = [](Point3D){ return 0.; };
    Func fNeu   = [](const Point3D & p) { return M_PI * std::sin(M_PI * p.y()); };
    Func SS     = [](const Point3D & p) { return 2 * M_PI * M_PI * std::sin(M_PI * p.x()) * std::sin(M_PI * p.y()); };
    Func u_ex   = [](const Point3D & p) { return std::sin(M_PI * p.x()) * std::sin(M_PI * p.y()); };

//    Func fNeu   = [](const Point3D & p) { return - p.y() - 1.; };
//    Func SS     = [](Point3D) { return 0.; };
//    Func u_ex   = [](const Point3D & p) { return p.x() * p.y() + p.x() - p.y(); };

    BoundaryConditions::BorderBC leftBC (BorderLabel::Left, Dirichlet, u_ex );
//    BoundaryConditions::BorderBC rightBC(BorderLabel::Right, Neumann, fNeu );
    BoundaryConditions::BorderBC rightBC(BorderLabel::Right, Dirichlet, u_ex );
    BoundaryConditions::BorderBC backBC (BorderLabel::Back, Dirichlet, u_ex );
    BoundaryConditions::BorderBC frontBC(BorderLabel::Front, Dirichlet, u_ex );
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

    Rigid_Mesh myrmesh(mesh, propMap);

    Pb * darcy(nullptr);
    darcy = new DarcyPb(dataPtr->getSolverType(), myrmesh, BC, SS, dataPtr);

    darcy->assemble();
    darcy->solve();

    Vector u_h_ex = evaluateMatrix(myrmesh, u_ex);
    Vector u_h = darcy->getSolver().getSolution();
    Vector u_diff = u_h_ex - u_h ;

    Quadrature q(myrmesh);

    const Real err = q.L2NormMatrix( u_h_ex - darcy->getSolver().getSolution() );
    std::cout << std::setprecision(15) << "L2 norm: " << err << std::endl;

#ifdef FVCODE3D_EXPORT
    ExporterVTU exporter;
    exporter.exportSolution(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_solution.vtu", darcy->getSolver().getSolution());
    exporter.exportSolution(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_solution_diff.vtu",
    u_diff.cwiseAbs());
    exporter.exportWithProperties(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_APP.vtu", Aperture | Permeability | Porosity);
#endif // FVCODE3D_EXPORT

    delete darcy;
    delete importer;

    if(std::fabs(err - 0.000777859293974688) <= 1e-13)
    {
        return EXIT_SUCCESS;
    }
    else
    {
        return EXIT_FAILURE;
    }
}
