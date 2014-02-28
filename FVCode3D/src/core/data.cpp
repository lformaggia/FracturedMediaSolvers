/*!
 *	@file chrono.cpp
 *	@brief This class implements a simple chronometer (definitions).
 */

#include "core/data.hpp"

Data::Data(const std::string dataFileName)
{
	GetPot dataFile((dataFileName).c_str());

	M_meshDir = dataFile("mesh/mesh_dir", "./data/");
	M_meshFile = dataFile("mesh/mesh_file", "grid.grid");

	M_outputDir = dataFile("output/output_dir", "./result/");
	M_outputFile = dataFile("output/output_file", "sol");

	M_type = static_cast<ProblemType>( dataFile("problem/type", steady) );

	M_mobility = dataFile("fluid/mobility", 1.);

	M_theta = dataFile("bc/theta", 0.);

	M_verbose = dataFile("miscellaneous/verbose", true);
}

void Data::showMe( std::ostream & output ) const
{
	output << "-----------------------------" << std::endl;

	output << "Mesh Directory: " << M_meshDir << std::endl;
	output << "Mesh Filename: " << M_meshFile << std::endl;

	output << "Output Directory: " << M_outputDir << std::endl;
	output << "Output Filename: " << M_outputFile << std::endl;

	output << "Type of the problem: ";
	switch(M_type)
	{
		case steady:
			output << "steady" << std::endl;
			break;
		case unsteady:
			output << "unsteady" << std::endl;
			break;
		default:
			output << "unsteady" << std::endl;
			break;
	}

	output << "Mobility: " << M_mobility << std::endl;

	output << "Theta: " << M_theta << std::endl;

	output << "Verbose: " << M_verbose << std::endl;

	output << "-----------------------------" << std::endl;
}
