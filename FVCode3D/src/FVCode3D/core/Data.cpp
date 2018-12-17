/*!
 *  @file data.cpp
 *  @brief This class handles the data for the solver (definitions).
 */

#include <FVCode3D/core/Data.hpp>
#include <FVCode3D/solver/SolverHandler.hpp>
#include <algorithm>
#include <tuple>
#include <vector>
#include <string>
#include <map>

namespace FVCode3D
{

Data::Data():
    M_meshDir("./data/"), M_meshFile("grid.fvg"), M_meshExt(".fvg"), M_meshType(forSolver),
    M_outputDir("./results/"), M_outputFile("sol"),
    M_Lx(2.), M_Ly(1.), M_Lz(1.), M_Nx(10), M_Ny(5), M_Nz(5),
    M_Sx(0.), M_Sy(0.), M_Sz(0.),
    M_Rz(0.),
    M_noise(false),
    M_noiseOn(NoiseOn::Matrix),
    M_meanNDist(0.), M_stDevNDist(1.),
    M_numet(FV), M_lumpedMim(false),
    M_problemType(steady), M_fracturesOn(true), M_ssOn(SourceSinkOn::None),
    M_setFracturesPressure(false), M_fracturesPressure(0.),
    M_nbSubRegions(1),
    M_nbTimeStepSteadyState(0), M_tolSteadyState(1e-8),
    M_permeabilityType("ScalarPermeability"),
    M_permMatrix(0.), M_poroMatrix(0.),
    M_permFrac(0.), M_poroFrac(0.), M_aperFrac(0.),
    M_initTime(0.), M_endTime(1.), M_timeStep(0.1),
    M_mobility(1.), M_compressibility(1.),
    M_SolverPolicy(SolverPolicy::Direct),
    M_maxIt(1.), M_tol(0.),M_restart(100u),
    M_theta(0.), M_verbose(true)
{}

Data::Data(const std::string dataFileName):Data()
{
  this->load(dataFileName);
}

void Data::load(const std::string dataFileName)
{
    EnumParser<MeshFormatType> parserMeshType;
    EnumParser<NumericalMethodType> parserNumericalMethodType;
    EnumParser<SolverPolicy> parserSolverPolicy;
    EnumParser<ProblemType> parserProblemType;
    EnumParser<NoiseOn> parserNoiseOn;
    EnumParser<SourceSinkOn> parserSourceSinkOn;

    std::ifstream file(dataFileName.c_str());
    if(file.good())
    {
        file.close();
    }
    else
    {
        file.close();
        throw std::runtime_error("Error: data file not opened.");
    }

    GetPot dataFile((dataFileName).c_str());

    M_meshDir = dataFile("mesh/mesh_dir", "./data/");
    M_meshFile = dataFile("mesh/mesh_file", "grid1.grid");
    M_meshExt = dataFile("mesh/mesh_type", ".grid");
    M_meshType = parserMeshType.parse( M_meshExt );

    M_outputDir = dataFile("output/output_dir", "./results/");
    M_outputFile = dataFile("output/output_file", "sol");

    M_Lx = dataFile("domain/Lx", 2.);
    M_Ly = dataFile("domain/Ly", 1.);
    M_Lz = dataFile("domain/Lz", 1.);
    M_Nx = dataFile("domain/Nx", 10);
    M_Ny = dataFile("domain/Ny", 5);
    M_Nz = dataFile("domain/Nz", 5);
    M_Sx = dataFile("domain/Sx", 0.);
    M_Sy = dataFile("domain/Sy", 0.);
    M_Sz = dataFile("domain/Sz", 0.);
    M_Rz = dataFile("domain/Rz", 0.);
    M_noise = static_cast<bool>(dataFile("domain/noise", 0));
    M_noiseOn = parserNoiseOn.parse( dataFile("domain/noiseOn", "matrix") );
    M_meanNDist = dataFile("domain/meanN", 0.);
    M_stDevNDist = dataFile("domain/stDevN", 1.);

    M_numet = parserNumericalMethodType.parse( dataFile("numet/method", "FV") );
    M_SolverPolicy = parserSolverPolicy.parse( dataFile("solver/policy", "Direct") );
    M_lumpedMim = static_cast<bool>(dataFile("mimetic/lumped", 0));

    M_problemType = parserProblemType.parse( dataFile("problem/type", "steady") );
    M_fracturesOn = static_cast<bool>(dataFile("problem/fracturesOn", 1));
    M_ssOn = parserSourceSinkOn.parse( dataFile("problem/sourceOn", "none") );

    M_setFracturesPressure = static_cast<bool>(dataFile("problem/fracPressOn", 0));
    M_fracturesPressure = dataFile("problem/fracPress", 0.);

    M_nbSubRegions = dataFile("msr/nbSubRegions", 1 );
    M_nbTimeStepSteadyState = dataFile("msr/nbStep", 0 );
    M_tolSteadyState = dataFile("msr/tol", 1e-8 );

    M_permeabilityType = dataFile("problem/perm_size", "ScalarPermeability");
    M_permMatrix = dataFile("problem/perm_matrix", 1e2);
    M_poroMatrix = dataFile("problem/poro_matrix", 0.25);
    M_permFrac = dataFile("problem/perm_frac", 1e5);
    M_poroFrac = dataFile("problem/poro_frac", 1.);
    M_aperFrac = dataFile("problem/aper_frac", 0.1);

    M_initTime = dataFile("problem/initial_time", 0.);
    M_endTime = dataFile("problem/end_time", 1.);
    M_timeStep = dataFile("problem/time_step", 0.1);

    M_mobility = dataFile("fluid/mobility", 1.);
    M_compressibility = dataFile("fluid/compressibility", 1.);

    M_solverType = dataFile("solver/type", "EigenUmfPack");
    M_maxIt = dataFile("solver/iterative/maxIt", 1000);
    M_tol = dataFile("solver/iterative/tolerance", 1e-4);
    M_precon = dataFile("solver/iterative/preconditioner", "ILU");
    M_lumpedMim=static_cast<bool>(dataFile("solver/iterative/lumped",0));
    M_restart=dataFile("solver/iterative/restart",100);

    M_theta = dataFile("bc/theta", 0.);

    M_verbose = static_cast<bool>(dataFile("miscellaneous/verbose", 1));

    std::tie(this->M_bulkPermeabilityData,this->M_fracturePermeabilityData)
    =Utility::readPermeabilityData(dataFile);
    this->M_bulkData=Utility::readBulkData(dataFile);
    this->M_fractureData=Utility::readFractureData(dataFile);
    this->M_dumpMatrix=static_cast<bool>(dataFile("miscellaneous/dumpMatrices", 0));
}

void Data::setMeshExtension(const std::string ext)
{
    EnumParser<MeshFormatType> parserMeshType;
    M_meshExt = ext;
    M_meshType = parserMeshType.parse( M_meshExt );
}

void Data::setMeshType(const MeshFormatType type)
{
    switch(M_meshType)
    {
        case TPFA:
            M_meshExt = ".grid";
            break;
        case forSolver:
            M_meshExt = ".fvg";
            break;
        case Medit:
            M_meshExt = ".mesh";
            break;
        case TetGen:
            M_meshExt = ".node";
            break;
        case OpenFOAM:
            M_meshExt = ".foam";
            break;
        default:
            throw std::runtime_error("Error: the mesh type set does not exist.");
            break;
    }

    M_meshType = type;
}

void Data::showMe( std::ostream & output ) const
{
    output << "-----------------------------" << std::endl;

    output << "Mesh Directory: " << M_meshDir << std::endl;
    output << "Mesh Filename: " << M_meshFile << std::endl;
    output << "Mesh Extension: " << M_meshExt << std::endl;
    output << "Mesh Type: ";
    switch(M_meshType)
    {
        case TPFA:
            output << "TPFA" << std::endl;
            break;
        case forSolver:
            output << "forSolver" << std::endl;
            break;
        case Medit:
            output << "Medit" << std::endl;
            break;
        case TetGen:
            output << "TetGen" << std::endl;
            break;
        case OpenFOAM:
            output << "OpenFOAM" << std::endl;
            break;
        default:
            exit(0);
            break;
    }

    output << "Output Directory: " << M_outputDir << std::endl;
    output << "Output Filename: " << M_outputFile << std::endl;

    output << "Numerical Method: ";
    switch(M_numet)
    {
        case FV:
            output << "Finite Volume" << std::endl;
            break;
        case MFD:
            output << "Mimetic Finite Difference" << std::endl;
            break;
        default:
            exit(0);
            break;
    }

    if(M_numet == MFD)
    {
        output << "Lumped: " << M_lumpedMim << std::endl;
    }

    output << "Type of the problem: ";
    switch(M_problemType)
    {
        case steady:
            output << "steady" << std::endl;
            break;
        case pseudoSteady:
            output << "pseudo-steady" << std::endl;
            break;
        default:
            exit(0);
            break;
    }
    output << "Fractures On: " << M_fracturesOn << std::endl;

    output << "Source/Sink on: ";
    switch(M_ssOn)
    {
        case SourceSinkOn::Matrix:
            output << "matrix" << std::endl;
            break;
        case SourceSinkOn::Fractures:
            output << "fractures" << std::endl;
            break;
        case SourceSinkOn::Both:
            output << "all" << std::endl;
            break;
        case SourceSinkOn::None:
            output << "none" << std::endl;
            break;
        default:
            exit(0);
            break;
    }

    output << "Fix pressure in fractures: " << M_setFracturesPressure << std::endl;

    if(M_setFracturesPressure)
    {
        output << "Fractures pressure: " << M_fracturesPressure << std::endl;
    }

    output << "Size of the permeability: " << M_permeabilityType << std::endl;

    if(M_problemType == pseudoSteady)
    {
        output << "Initial time: " << M_initTime << std::endl;
        output << "End time: " << M_endTime << std::endl;
        output << "Time step: " << M_timeStep << std::endl;

        output << "Compressibility: " << M_compressibility << std::endl;
    }

    output << "Mobility: " << M_mobility << std::endl;

    output << "Solver: " << M_solverType << std::endl;
    {
      auto SolverPtr = SolverHandler::Instance().getProduct(M_solverType);
      bool testIterative = dynamic_cast<IterativeSolver*>(SolverPtr.get());
      if(testIterative)
        {
          output << "# max iter: " << M_maxIt << std::endl;
          output << "Tolerance: " << M_tol << std::endl;
          if(this->M_lumpedMim)
            output<<"Using lumped inner product matrix for approximating Shur complement"<<std::endl;
          else
            output<<"Using diagonal of inner product matrix for approximating Shur complement"<<std::endl;
          bool testGMRES = dynamic_cast<imlGMRES*>(SolverPtr.get()) || dynamic_cast<imlFGMRES*>(SolverPtr.get());
          if (testGMRES)
            output<<"Restart level:  "<<this->M_restart<<std::endl;
          output << "Preconditioner: " << M_precon << std::endl;
        }
    }
    output << "Theta: " << M_theta << std::endl;

    output << "Verbose: " << M_verbose << std::endl;

    output << "-----------------------------" << std::endl;

    std::map<int,std::string> permName{{0,"Scalar"},{1,"Diagonal"},{2,"SymTensor"}};

    output<<"  MATRIX PROPERTIES\n";
    for (auto i=0u; i< this->M_bulkPermeabilityData.size(); ++i)
      {
        auto & g = this->M_bulkPermeabilityData[i];
        output<<"Zone "<<std::get<0>(g)<<" ";
        output<<"Permeability type "<<permName[std::get<0>(std::get<1>(g))]<<" Values: ";
        for (auto & p : std::get<1>(std::get<1>(g))) output<<p<<" ";
        output<<"\n";
        output<<"Porosity: "<<std::get<1>(this->M_bulkData[i])<<"\n";
      }
    output<<"FRACTURE PROPERTIES\n";
    for (auto i=0u; i< this->M_fracturePermeabilityData.size(); ++i)
       {
         auto & g = this->M_fracturePermeabilityData[i];
         output<<"Zone "<<std::get<0>(g)<<" ";
         output<<"Permeability type "<<permName[std::get<0>(std::get<1>(g))]<<" Values: ";
         for (auto & p : std::get<1>(std::get<1>(g))) output<<p<<" ";
         output<<"\n";
         output<<"Porosity: "<<std::get<1>(this->M_fractureData[i])<<"Aperture: "<<std::get<2>(this->M_fractureData[i])<<"\n";
       }

}

std::tuple<double,double> Data::fracturePorosityAndAperture(int const & zoneNumber) const
{
    return Utility::getFractureData(this->M_fractureData,zoneNumber);
}
  
Real Data::getFractureAperture(int zoneNumber) const
{
    auto r = Data::fracturePorosityAndAperture(zoneNumber);
    return std::get<1>(r);
}

Real Data::getFracturePorosity(int zoneNumber) const
{
    auto r = Data::fracturePorosityAndAperture(zoneNumber);
    return std::get<0>(r);
}

Real Data::getMatrixPorosity(int zoneNumber) const
{
    return Utility::getBulkData(this->M_bulkData,zoneNumber);
}

Utility::Permeability
Data::getFracturePermeability(int zoneNumber) const
{
    return Utility::getPermeability(this->M_fracturePermeabilityData,zoneNumber);
}

Utility::Permeability
Data::getMatrixPermeability(int zoneNumber) const
{
    return Utility::getPermeability(this->M_bulkPermeabilityData,zoneNumber);
}

template<class T>
EnumParser<T>::EnumParser() = default;

template<>
EnumParser<Data::NumericalMethodType>::EnumParser()
{
    M_enumMap[ toUpper( "FV" ) ] = Data::NumericalMethodType::FV;
    M_enumMap[ toUpper( "MFD" ) ] = Data::NumericalMethodType::MFD;
}

template<>
EnumParser<Data::SolverPolicy>::EnumParser()
{
    M_enumMap[ toUpper( "Direct" ) ] = Data::SolverPolicy::Direct;
    M_enumMap[ toUpper( "Iterative" ) ] = Data::SolverPolicy::Iterative;
}

template<>
EnumParser<Data::ProblemType>::EnumParser()
{
    M_enumMap[ toUpper( "steady" ) ] = Data::ProblemType::steady;
    M_enumMap[ toUpper( "pseudoSteady" ) ] = Data::ProblemType::pseudoSteady;
}

template<>
EnumParser<Data::MeshFormatType>::EnumParser()
{
    M_enumMap[ toUpper( ".grid" ) ] = Data::MeshFormatType::TPFA;
    M_enumMap[ toUpper( ".fvg" ) ] = Data::MeshFormatType::forSolver;
    M_enumMap[ toUpper( ".mesh" ) ] = Data::MeshFormatType::Medit;
    M_enumMap[ toUpper( ".node" ) ] = Data::MeshFormatType::TetGen;
    M_enumMap[ toUpper( ".foam" ) ] = Data::MeshFormatType::OpenFOAM;
}

template<>
EnumParser<Data::NoiseOn>::EnumParser()
{
    M_enumMap[ toUpper( "matrix" ) ] = Data::NoiseOn::Matrix;
    M_enumMap[ toUpper( "fractures" ) ] = Data::NoiseOn::Fractures;
    M_enumMap[ toUpper( "all" ) ] = Data::NoiseOn::All;
}

template<>
EnumParser<Data::SourceSinkOn>::EnumParser()
{
    M_enumMap[ toUpper( "matrix" ) ] = Data::SourceSinkOn::Matrix;
    M_enumMap[ toUpper( "fractures" ) ] = Data::SourceSinkOn::Fractures;
    M_enumMap[ toUpper( "all" ) ] = Data::SourceSinkOn::Both;
    M_enumMap[ toUpper( "none" ) ] = Data::SourceSinkOn::None;
}

template<>
EnumParser<Data::PermeabilityType>::EnumParser()
{
    M_enumMap[ toUpper( "scalar" ) ] = Data::PermeabilityType::Scalar;
    M_enumMap[ toUpper( "diagonal" ) ] = Data::PermeabilityType::Diagonal;
    M_enumMap[ toUpper( "symtensor" ) ] = Data::PermeabilityType::SymTensor;
    M_enumMap[ toUpper( "arbitrarytensor" ) ] = Data::PermeabilityType::ArbitraryTensor;
}
// Definition of the global namespace variable to hold general parameters for the code
DataPtr_Type dataPtr(new Data());
} // namespace FVCode3D
