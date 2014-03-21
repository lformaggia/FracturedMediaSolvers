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
	 * This enumerator allows to select the format type of the mesh: TPFA(.grid) or forSolver(.fvg)
	 */
	enum MeshFormatType
	{
		TPFA			= 0,
		forSolver		= 1
	};

	//! @name Constructors
	//@{

	//! Default constructor
	/*!
	 * Initialize the data
	 */
	Data();

	//! Constructor from file
	/*!
	 * @param dataFileName name of the data file
	 */
	Data(const std::string dataFileName);

	//! Copy constructor
	/*!
	 * @param data reference of a Data
	 */
	Data(const Data & data);

	//@}

	//! @name Get methods
	//@{

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

	//! Get the extension of the mesh format type
	/*!
	 * @return the extension of the format type of the mesh
	 */
	const std::string getMeshExtension() const { return M_meshExt; }

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

	//! Get the initial time
	/*!
	 * @return the initial time
	 */
	Real getInitialTime() const { return M_initTime; }

	//! Get the end time
	/*!
	 * @return the end time
	 */
	Real getEndTime() const { return M_endTime; }

	//! Get the time step
	/*!
	 * @return the time step
	 */
	Real getTimeStep() const { return M_timeStep; }

	//! Get mobility of the fluid
	/*!
	 * @return the mobility of the fluid
	 */
	Real getMobility() const { return M_mobility; }

	//! Get compressibility
	/*!
	 * @return the compressibility
	 */
	Real getCompressibility() const { return M_compressibility; }

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

	//@}

	//! @name Set methods
	//@{

	//! Set the directory of the mesh
	/*!
	 * @param dir the string that contains the directory of the mesh
	 */
	void setMeshDir(const std::string dir) { M_meshDir = dir; }

	//! Set the file of the mesh
	/*!
	 * @param file the string that contains the file of the mesh
	 */
	void setMeshFile(const std::string file) { M_meshFile = file; }

	//! Set the extension of the mesh format type
	/*!
	 * @param ext the extension of the format type of the mesh
	 * @post modifies the mesh type
	 */
	void setMeshExtension(const std::string ext);

	//! Set the format type of the mesh
	/*!
	 * @param type the format type of the mesh
	 * @post modifies the mesh extension
	 */
	void setMeshType(const MeshFormatType type);
	
	//! Set the output directory
	/*!
	 * @param dir the string that contains the output directory
	 */
	void setOutputDir(const std::string dir) { M_outputDir = dir; }

	//! Set the output file name
	/*!
	 * @param file the string that contains the output file name
	 */
	void setOutputFile(const std::string file) { M_outputFile = file; }

	//! Set the problem type
	/*!
	 * @param type the type of the problem
	 */
	void setProblemType(const ProblemType type) { M_problemType = type; };

	//! Enable or disable the fractures
	/*!
	 * @param fracture true to enable the fractures
	 */
	void fractureOn(bool fracture) { M_fracturesOn = fracture; }

	//! Set the initial time
	/*!
	 * @param time the initial time
	 */
	void setInitialTime(const Real time) { M_initTime = time; }

	//! Set the final time
	/*!
	 * @param time the final time
	 */
	void setEndTime(const Real time) { M_endTime = time; }

	//! Set the time step
	/*!
	 * @param time the time step
	 */
	void setTimeStep(const Real time) { M_timeStep = time; }

	//! Set mobility of the fluid
	/*!
	 * @param mobility the mobility of the fluid
	 */
	void setMobility(const Real mobility) { M_mobility = mobility; }

	//! Set compressibility
	/*!
	 * @param compressibility the compressibility
	 */
	void setCompressibility(const Real compressibility) { M_compressibility = compressibility; }

	//! Set theta angle
	/*!
	 * @param theta the angle theta used to identify the BC ids
	 */
	void setTheta(const Real theta) { M_theta = theta; }

	//! Enable or disable verbosity
	/*!
	 * @param verbose true to enable verbosity
	 */
	void verbose(const bool verbose) { M_verbose = verbose; }

	//@}

	//! Show a summary of the data
	/*!
	 * @param output output stream. Default = std::cout
	 */
	void showMe( std::ostream & output = std::cout ) const;

	//! Destructor
	~Data() {}

protected:

	//! Mesh directory
	std::string M_meshDir;
	//! Mesh file
	std::string M_meshFile;
	//! Mesh file
	std::string M_meshExt;
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
	//! Initial time, used for time-dependent problem
	Real M_initTime;
	//! Final time, used for time-dependent problem
	Real M_endTime;
	//! Time step, used for time-dependent problem
	Real M_timeStep;
	//! Mobility of the fluid
	Real M_mobility;
	//! Compressibility
	Real M_compressibility;
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
