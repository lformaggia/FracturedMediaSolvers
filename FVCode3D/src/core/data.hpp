/*
 *	@file data.hpp
 *	@brief This class handles the data for the solver.
 */

#ifndef DATA_HPP_
#define DATA_HPP_

#include "core/TypeDefinition.hpp"

//! Class that handles the information needed for the solver
/*!
 * @class Data
 * This class takes care of the information needed for the solver, as boundary conditions and the source/sink term.
 * It also handles the name of the input/output directory and the name of the input/output file,
 * the shared library name and miscellaneous options.
 */
class Data{
public:

	//! Define the problem type
	/*!
	 * @enum ProblemType
	 * This enumerator allows to select the problem type: steady or pseudo-steady state
	 */
	enum ProblemType
	{
		steady			= 0,
		pseudoSteady	= 1
	};

	//! Define the format of the input mesh
	/*!
	 * @enum MeshFormatType
	 * This enumerator allows to select the format type of the mesh: TPFA(old) or forSolver(new)
	 */
	enum MeshFormatType
	{
		TPFA			= 0,
		forSolver		= 1
	};

	//! Constructor
	/*!
	 * @param dataFileName name of the data file
	 */
	Data(const std::string dataFileName);

	//! Get the directory of the mesh
	/*!
	 * @return the string that contains the directory of the mesh
	 */
	const std::string getMeshDir() const { return M_meshDir; }

	//! Get the file of the mesh
	/*!
	 * @return the string that contains the file of the mesh
	 */
	const std::string getMeshFile() const { return M_meshFile; }

	//! Get the format type of the mesh
	/*!
	 * @return the format type of the mesh
	 */
	MeshFormatType getMeshType() const { return M_meshType; }
	
	//! Get the output directory
	/*!
	 * @return the string that contains the output directory
	 */
	const std::string getOutputDir() const { return M_outputDir; }

	//! Get the output file name
	/*!
	 * @return the string that contains the output file name
	 */
	const std::string getOutputFile() const { return M_outputFile; }

	//! Get the problem type
	/*!
	 * @return the type of the problem
	 */
	ProblemType getProblemType() const { return M_problemType; }

	//! Test if fractures are enabled
	/*!
	 * @return true if the fractures are enabled
	 */
	bool fractureOn() const { return M_fracturesOn; }

	//! Get mobility of the fluid
	/*!
	 * @return the mobility of the fluid
	 */
	Real getMobility() const { return M_mobility; }

	//! Get theta angle
	/*!
	 * @return the angle theta used to identify the BC ids
	 */
	Real getTheta() const { return M_theta; }

	//! Verbose
	/*!
	 * @return true if verbose
	 */
	bool verbose() const { return M_verbose; }

	//! Show a summary of the data
	/*!
	 * @param output output stream. Default = std::cout
	 */
	void showMe( std::ostream & output = std::cout ) const;

	//! Destructor
	~Data() {}

protected:

	//! No empty constructor
	Data();

	//! No copy constructor
	Data(const Data &);

	//! Mesh directory
	std::string M_meshDir;
	//! Mesh file
	std::string M_meshFile;
	//! Format Type
	MeshFormatType M_meshType;
	//! Output directory
	std::string M_outputDir;
	//! Output file
	std::string M_outputFile;
	//! Type of the problem
	ProblemType M_problemType;
	//! Enable or disable fractures
	bool M_fracturesOn;
	//! Mobility of the fluid
	Real M_mobility;
	//! Angle used to rotate along z-axis the domain. It is used only to compute the normal for detecting BC!
	Real M_theta;
	//! Verbose
	bool M_verbose;

};

template <class T>
class EnumParser
{
public:
	EnumParser(){};

	T parse(const std::string & str) const
	{ 
		typename std::map<std::string, T>::const_iterator it = enumMap.find(str);
		if (it == enumMap.end())
		{	
			std::cerr << "This enum does not exist!" << std::endl;
			exit(0);
		}
		return it->second;
	}
private:
	std::map<std::string, T> enumMap;
};

template<>
EnumParser<Data::ProblemType>::EnumParser();

template<>
EnumParser<Data::MeshFormatType>::EnumParser();

#endif /* DATA_HPP_ */