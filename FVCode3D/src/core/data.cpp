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
	M_meshFile = dataFile("mesh/mesh_file", "grid.grid");
	M_meshType = parserMeshType.parse( dataFile("mesh/mesh_type", "forSolver") );

	M_outputDir = dataFile("output/output_dir", "./result/");
	M_outputFile = dataFile("output/output_file", "sol");

	M_problemType = parserProblemType.parse( dataFile("problem/type", "steady") );

	M_fracturesOn = dataFile("problem/fracturesOn", true);

	M_mobility = dataFile("fluid/mobility", 1.);

	M_theta = dataFile("bc/theta", 0.);

	M_verbose = dataFile("miscellaneous/verbose", true);
}

void Data::showMe( std::ostream & output ) const
{
	output << "-----------------------------" << std::endl;

	output << "Mesh Directory: " << M_meshDir << std::endl;
	output << "Mesh Filename: " << M_meshFile << std::endl;
	output << "Mesh Type: ";
	switch(M_meshType)
	{
		case steady:
			output << "TPFA" << std::endl;
			break;
		case pseudoSteady:
			output << "forSolver" << std::endl;
			break;
		default:
			output << "forSolver" << std::endl;
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
			output << "pseudo-steady" << std::endl;
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
    enumMap["TPFA"] = Data::MeshFormatType::TPFA;
    enumMap["forSolver"] = Data::MeshFormatType::forSolver;
}