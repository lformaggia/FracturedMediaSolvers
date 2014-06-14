/*!
 *	@file chrono.cpp
 *	@brief This class implements a simple chronometer (definitions).
 */

#include "core/data.hpp"

namespace FVCode3D
{

Data::Data():
	M_meshDir("./data/"), M_meshFile("grid.fvg"), M_meshExt(".fvg"), M_meshType(forSolver),
	M_outputDir("./results/"), M_outputFile("sol"),
	M_Lx(2.), M_Ly(1.), M_Lz(1.), M_Nx(10), M_Ny(5), M_Nz(5),
	M_problemType(steady), M_fracturesOn(true), M_ssOn(None),
	M_setFracturesPressure(false), M_fracturesPressure(0.),
	M_permMatrix(0.), M_poroMatrix(0.),
	M_permFrac(0.), M_poroFrac(0.), M_aperFrac(0.),
	M_initTime(0.), M_endTime(1.), M_timeStep(0.1),
	M_mobility(1.), M_compressibility(0.), M_theta(0.), M_verbose(true)
{}

Data::Data(const std::string dataFileName)
{
	EnumParser<MeshFormatType> parserMeshType;
	EnumParser<ProblemType> parserProblemType;
	EnumParser<SourceSinkOn> parserSourceSinkOn;
  
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

	M_problemType = parserProblemType.parse( dataFile("problem/type", "steady") );
	M_fracturesOn = static_cast<bool>(dataFile("problem/fracturesOn", 1));
	M_ssOn = parserSourceSinkOn.parse( dataFile("problem/sourceOn", "none") );

	M_setFracturesPressure = static_cast<bool>(dataFile("problem/fracPressOn", 0));;
	M_fracturesPressure = dataFile("problem/fracPress", 0.);

	M_permMatrix = dataFile("problem/perm_matrix", 1e2);
	M_poroMatrix = dataFile("problem/poro_matrix", 0.25);
	M_permFrac = dataFile("problem/perm_frac", 1e5);
	M_poroFrac = dataFile("problem/poro_frac", 1.);
	M_aperFrac = dataFile("problem/aper_frac", 0.1);

	M_initTime = dataFile("problem/initial_time", 0.);
	M_endTime = dataFile("problem/end_time", 1.);
	M_timeStep = dataFile("problem/time_step", 0.1);

	M_mobility = dataFile("fluid/mobility", 1.);
	M_compressibility = dataFile("fluid/compressibility", 0.);

	M_theta = dataFile("bc/theta", 0.);

	M_verbose = static_cast<bool>(dataFile("miscellaneous/verbose", 1));
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
		default:
			exit(0);
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
		default:
			exit(0);
			break;
	}

	output << "Output Directory: " << M_outputDir << std::endl;
	output << "Output Filename: " << M_outputFile << std::endl;

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
		case Matrix:
			output << "matrix" << std::endl;
			break;
		case Fractures:
			output << "fractures" << std::endl;
			break;
		case Both:
			output << "all" << std::endl;
			break;
		case None:
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

	if(M_problemType == pseudoSteady)
	{
		output << "Initial time: " << M_initTime << std::endl;
		output << "End time: " << M_endTime << std::endl;
		output << "Time step: " << M_timeStep << std::endl;
	}

	output << "Mobility: " << M_mobility << std::endl;

	if(M_problemType == pseudoSteady)
		output << "Compressibility: " << M_compressibility << std::endl;

	output << "Theta: " << M_theta << std::endl;

	output << "Verbose: " << M_verbose << std::endl;

	output << "-----------------------------" << std::endl;
}

template<>
EnumParser<Data::ProblemType>::EnumParser()
{
    enumMap["steady"] = Data::ProblemType::steady;
    enumMap["pseudoSteady"] = Data::ProblemType::pseudoSteady;
}

template<>
EnumParser<Data::MeshFormatType>::EnumParser()
{
    enumMap[".grid"] = Data::MeshFormatType::TPFA;
    enumMap[".fvg"] = Data::MeshFormatType::forSolver;
    enumMap[".mesh"] = Data::MeshFormatType::Medit;
}

template<>
EnumParser<Data::SourceSinkOn>::EnumParser()
{
    enumMap["matrix"] = Data::SourceSinkOn::Matrix;
    enumMap["fractures"] = Data::SourceSinkOn::Fractures;
    enumMap["all"] = Data::SourceSinkOn::Both;
    enumMap["none"] = Data::SourceSinkOn::None;
}

} // namespace FVCode3D
