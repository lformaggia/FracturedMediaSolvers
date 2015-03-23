//#define NDEBUG // add this macro to disable asserts
#include <cassert>
#include <iostream>
#include <chrono>
#include <vector>
#include <fstream>

#include <FVCode3D/core/TypeDefinition.hpp>
#include <FVCode3D/core/Data.hpp>
#include <FVCode3D/mesh/Mesh3D.hpp>
#include <FVCode3D/property/Properties.hpp>
#include <FVCode3D/property/Permeability.hpp>
#include <FVCode3D/import/Import.hpp>
#include <FVCode3D/export/ExportVTU.hpp>
#include <FVCode3D/utility/Converter.hpp>

using namespace FVCode3D;

int main(int argc, char * argv[])
{
    GetPot command_line(argc,argv);
    const std::string dataFileName = command_line.follow("data.txt", 2, "-f", "--file");

    std::cout << "Read Data..." << std::flush;
    DataPtr_Type dataPtr(new Data(dataFileName));
    dataPtr->verbose(true);
    std::cout << " done." << std::endl;

    std::cout << std::endl;
    dataPtr->showMe();
    std::cout << std::endl;


    std::cout << "Define Mesh and Properties..." << std::flush;
    Mesh3D mesh;
    PropertiesMap propMap(dataPtr->getMobility());
    std::cout << " done." << std::endl;


    std::cout << "Create Importer..." << std::flush;
    Importer * importer = 0;
    if(dataPtr->getMeshType() == Data::MeshFormatType::Medit)
        importer = new ImporterMedit(dataPtr->getMeshDir() + dataPtr->getMeshFile(), mesh, propMap);
    else if(dataPtr->getMeshType() == Data::MeshFormatType::TetGen)
        importer = new ImporterTetGen(dataPtr->getMeshDir() + dataPtr->getMeshFile(), mesh, propMap);
    else if(dataPtr->getMeshType() == Data::MeshFormatType::forSolver)
        importer = new ImporterForSolver(dataPtr->getMeshDir() + dataPtr->getMeshFile(), mesh, propMap);
    std::cout << " done." << std::endl;


    std::cout << "Import grid file..." << std::flush;
    importer->import(dataPtr->fractureOn());
    std::cout << " done." << std::endl << std::endl;


    std::cout << "Compute separated cells..." << std::flush;
    mesh.updateFacetsWithCells();
    std::cout << " done." << std::endl;


    std::cout << "Compute neighboring cells..." << std::flush;
    mesh.updateCellsWithNeighbors();
    std::cout << " done." << std::endl << std::endl;


    std::cout << "Set labels on boundary" << std::flush;
    if(dataPtr->getMeshType() == Data::MeshFormatType::Medit
        ||
       dataPtr->getMeshType() == Data::MeshFormatType::TetGen)
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


//    std::cout << "Generate Polyhedral grid..." << std::flush;

    // ...


    std::cout << "Set uniform properties..." << std::flush;
    std::shared_ptr<PermeabilityBase> matrixPerm( new PermeabilityScalar );
    std::shared_ptr<PermeabilityBase> fracturesPerm( new PermeabilityScalar );
    matrixPerm->setPermeability( dataPtr->getMatrixPermeability(), 0 );
    fracturesPerm->setPermeability( dataPtr->getFracturePermeability(), 0 );
    propMap.setPropertiesOnMatrix(mesh, dataPtr->getMatrixPorosity(), matrixPerm);
    propMap.setPropertiesOnFractures(mesh, dataPtr->getFractureAperture(), dataPtr->getFracturePorosity(), fracturesPerm);
    std::cout << " done." << std::endl << std::endl;

    std::cout << "Export..." << std::flush;
    saveAsOpenFOAMFormat(dataPtr->getOutputDir() + dataPtr->getOutputFile(), mesh);
    saveAsSolverFormat(dataPtr->getOutputDir() + dataPtr->getOutputFile() + ".fvg", mesh, propMap);
    std::cout << " done." << std::endl;

    delete importer;

    return 0;
}
