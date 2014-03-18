/*!
 *	@file chrono.cpp
 *	@brief This class implements a simple chronometer (definitions).
 */

#include "core/data.hpp"

Data::Data(const std::string dataFileName)
{
	EnumParser<MeshFormatType> parserMeshType;
	EnumParser<ProblemType> parserProblemType;
  
	GetPot dataFile((dataFileName).c_str());

	M_meshDir = dataFile("mesh/mesh_dir", "./data/");
	M_meshFile = dataFile("mesh/mesh_file", "grid1.grid");
	M_meshExt = dataFile("mesh/mesh_type", ".grid");
	M_meshType = parserMeshType.parse( M_meshExt );

	M_outputDir = dataFile("output/output_dir", "./results/");
	M_outputFile = dataFile("output/output_file", "sol");

	M_problemType = parserProblemType.parse( dataFile("problem/type", "steady") );

	M_fracturesOn = dataFile("problem/fracturesOn", true);

	M_mobility = dataFile("fluid/mobility", 1.);

	M_theta = dataFile("bc/theta", 0.);

	M_verbose = dataFile("miscellaneous/verbose", true);
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

	output << "Mobility: " << M_mobility << std::endl;

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
}
