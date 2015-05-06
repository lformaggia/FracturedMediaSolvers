#include <FVCode3D/FVCode3D.hpp>

//#define FVCODE3D_EXPORT

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

    Func fZero  = [](Point3D){ return 0.; };
    Func fOne   = [](Point3D){ return 1.; };
    Func SS     = [](Point3D){ return 0.; };

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

    Rigid_Mesh myrmesh(mesh, propMap);

    Pb * darcy(nullptr);
    darcy = new DarcyPb(dataPtr->getSolverType(), myrmesh, BC, SS, dataPtr);

    darcy->assemble();
    darcy->solve();

    Quadrature q(myrmesh);

    const Vector & uh = darcy->getSolver().getSolution();
    const Real norm_M = q.L2NormMatrix( uh );
    const Real norm_F = q.L2NormFractures( uh );
    const Real norm = std::sqrt(norm_M * norm_M + norm_F * norm_F);
    std::cout << std::setprecision(15) << "L2 norm: " << norm << std::endl;
    std::cout << std::setprecision(15) << "Matrix L2 norm: " << norm_M << std::endl;
    std::cout << std::setprecision(15) << "Fracture L2 norm: " << norm_F << std::endl;

#ifdef FVCODE3D_EXPORT
    ExporterVTU exporter;
    exporter.exportSolution(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_solution.vtu", darcy->getSolver().getSolution());
    exporter.exportSolutionOnFractures(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_solution_f.vtu", darcy->getSolver().getSolution());
    exporter.exportWithProperties(myrmesh, dataPtr->getOutputDir() + dataPtr->getOutputFile() + "_APP.vtu", Aperture | Permeability | Porosity);
#endif // FVCODE3D_EXPORT

    delete darcy;
    delete importer;

    if(std::fabs(norm - 1.68234751035572) <= 1e-13)
    {
        return EXIT_SUCCESS;
    }
    else
    {
        return EXIT_FAILURE;
    }
}
