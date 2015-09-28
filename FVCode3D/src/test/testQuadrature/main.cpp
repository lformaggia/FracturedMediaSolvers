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
    const Real kf = 1.e3;
    matrixPerm->setPermeability( 1., 0 );
    fracturesPerm->setPermeability( kf, 0 );
    constexpr Real aperture = 1.e-2;
    constexpr Real matrixPoro = 0.25;
    constexpr Real fracturesPoro = 1.;
    propMap.setPropertiesOnMatrix(mesh, matrixPoro, matrixPerm);
    propMap.setPropertiesOnFractures(mesh, aperture, fracturesPoro, fracturesPerm);

    Func fOne   = [](Point3D){ return 1.; };

    Rigid_Mesh myrmesh(mesh, propMap);

    Quadrature q(myrmesh);

    const Vector fM_h = q.cellIntegrateMatrix(fOne);
    const Vector fF_h = q.cellIntegrateFractures(fOne);
    const Real vM_h = fM_h.sum();
    const Real vF_h = fF_h.sum();
    constexpr Real vM = 4.;
    constexpr Real vF = 2. * aperture;
    const Real errvM = std::fabs(vM_h - vM) / vM;
    const Real errvF = std::fabs(vF_h - vF) / vF;
    std::cout << std::setprecision(15) << "Matrix volume: " << vM_h << std::endl;
    std::cout << std::setprecision(15) << "Fractures volume: " << vF_h << std::endl;
    std::cout << std::setprecision(15) << "Error matrix volume: " << errvM << std::endl;
    std::cout << std::setprecision(15) << "Error fracture volume: " << errvF << std::endl;
//    std::cout << std::setprecision(15) << "Total sum: " << vM_h + vF_h << std::endl;

    delete importer;

    if( errvM <= 1e-14 && errvF <= 1e-14 )
    {
        return EXIT_SUCCESS;
    }
    else
    {
        return EXIT_FAILURE;
    }
}
