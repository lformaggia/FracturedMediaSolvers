/*
 *  @file data.hpp
 *  @brief This class handles the data for the solver.
 */

#ifndef DATA_HPP_
#define DATA_HPP_

#include <FVCode3D/core/TypeDefinition.hpp>

namespace FVCode3D
{

//! Class that handles the information needed for the solver
/*!
 * @class Data
 * This class takes care of the information needed for the solver, as boundary conditions and the source/sink term.
 * It also handles the name of the input/output directory and the name of the input/output file,
 * the shared library name and miscellaneous options.
 */
class Data{
public:

    //! Define the numerical method type
    /*!
     * @enum NumericalMethodType
     * This enumerator allows to select the numerical method: Finite Volume or Mimetic Finite Difference
     */
    enum NumericalMethodType
    {
        FV          = 0,
        MFD         = 1
    };

    //! Define the problem type
    /*!
     * @enum ProblemType
     * This enumerator allows to select the problem type: steady or pseudo-steady state
     */
    enum ProblemType
    {
        steady          = 0,
        pseudoSteady    = 1
    };

    //! Define the format of the input mesh
    /*!
     * @enum MeshFormatType
     * This enumerator allows to select the format type of the mesh: TPFA(.grid) or forSolver(.fvg)
     */
    enum MeshFormatType
    {
        TPFA            = 0,
        forSolver       = 1,
        Medit           = 2
    };

    //! Define where to apply the source/sink term
    /*!
     * @enum SourceSinkOn
     * This enumerator allows to select where to apply the source/sink term: only on matrix, only on fractures, both or none.
     */
    enum SourceSinkOn
    {
        Matrix          = 0,
        Fractures       = 1,
        Both            = 2,
        None            = 3
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
    Data(const std::string dataFileName) throw();

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

    //! Get the domain dimension along x-axis
    /*!
     * @return the domain dimension along x-axis
     */
    Real getLx() const { return M_Lx; }

    //! Get the domain dimension along y-axis
    /*!
     * @return the domain dimension along y-axis
     */
    Real getLy() const { return M_Ly; }

    //! Get the domain dimension along z-axis
    /*!
     * @return the domain dimension along z-axis
     */
    Real getLz() const { return M_Lz; }

    //! Get the number of cells along x-axis
    /*!
     * @return the number of cells along x-axis
     */
    Real getNx() const { return M_Nx; }

    //! Get the number of cells along y-axis
    /*!
     * @return the number of cells along y-axis
     */
    Real getNy() const { return M_Ny; }

    //! Get the number of cells along z-axis
    /*!
     * @return the number of cells along z-axis
     */
    Real getNz() const { return M_Nz; }

    //! Get the shift domain along x-axis
    /*!
     * @return the shift domain along x-axis
     */
    Real getSx() const { return M_Sx; }

    //! Get the shift domain along y-axis
    /*!
     * @return the shift domain along y-axis
     */
    Real getSy() const { return M_Sy; }

    //! Get the shift domain z-axis
    /*!
     * @return the shift domain along z-axis
     */
    Real getSz() const { return M_Sz; }

    //! Get the numerical method used to solve the problem
    /*!
     * @return the numerical method used to solve the problem
     */
    NumericalMethodType getNumericalMethod() const { return M_numet; }

    //! Get if the stiffness matrix in the mimetic method is lumped
    /*!
     * @return lumped true if the stiffness matrix is lumped
     */
    bool getLumpedMimetic() const { return M_lumpedMim; }

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

    //! Get where the source/sink term is applied
    /*!
     * @return where the source/sink term is applied
     */
    SourceSinkOn getSourceSinkOn() const { return M_ssOn; }

    //! Test if the pressure value inside the fractures is fixed or not
    /*!
     * @return if true fix a pressure value inside the fractures
     */
    bool pressuresInFractures() const { return M_setFracturesPressure; }

    //! Get the pressure inside the fracture
    /*!
     * @return the pressure inside the fracture
     */
    Real getPressuresInFractures() const { return M_fracturesPressure; }

    //! Test if the multiple sub-regions method is activated or not
    /*!
     * @return if true the multiple sub-regions method is activated
     */
    bool MSROn() const { return M_MSR; }

    //! Get the number of sub-regions
    /*!
     * @return the number of sub-regions
     */
    UInt nbSubRegions() const { return M_nbSubRegions; }

    //! Get the number of time step to check for the steady state
    /*!
     * @return the number of time step to check for the steady state
     */
    UInt nbTimeStepSteadyState() const { return M_nbTimeStepSteadyState; }

    //! Get the tolerance to consider a time step at the steady state
    /*!
     * @return the tolerance to consider a time step at the steady state
     */
    Real tolSteadyState() const { return M_tolSteadyState; }

    //! Get the type of the permeability
    /*!
     * @return the type of the permeability
     */
    const std::string getPermeabilityType() const { return M_permeabilityType; }

    //! Get the permeability in the porous medium
    /*!
     * @return the permeability in the porous medium
     */
    Real getMatrixPermeability() const { return M_permMatrix; }

    //! Get the porosity in the porous medium
    /*!
     * @return the porosity in the porous medium
     */
    Real getMatrixPorosity() const { return M_poroMatrix; }

    //! Get the permeability in the fracture
    /*!
     * @return the permeability in the fracture
     */
    Real getFracturePermeability() const { return M_permFrac; }

    //! Get the porosity in the fracture
    /*!
     * @return the porosity in the fracture
     */
    Real getFracturePorosity() const { return M_poroFrac; }

    //! Get the aperture in the fracture
    /*!
     * @return the aperture in the fracture
     */
    Real getFractureAperture() const { return M_aperFrac; }

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

    //! Get solver type
    /*!
     * @return the solver type
     */
    const std::string getSolverType() const { return M_solverType; }

    //! Get maximum iterations of the iterative solver
    /*!
     * @return maximum iterations of the iterative solver
     */
    UInt getIterativeSolverMaxIter() const { return M_maxIt; }

    //! Get tolerance of the iterative solver
    /*!
     * @return tolerance of the iterative solver
     */
    Real getIterativeSolverTolerance() const { return M_tol; }

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
    void setMeshType(const MeshFormatType type) throw();

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

    //! Set the domain dimension along x-axis
    /*!
     * @param the domain dimension along x-axis
     */
    void setLx(const Real Lx) { M_Lx = Lx; }

    //! Set the domain dimension along y-axis
    /*!
     * @param the domain dimension along y-axis
     */
    void setLy(const Real Ly) { M_Ly = Ly; }

    //! Set the domain dimension along z-axis
    /*!
     * @param the domain dimension along z-axis
     */
    void setLz(const Real Lz) { M_Lz = Lz; }

    //! Set the number of cells along x-axis
    /*!
     * @param the number of cells along x-axis
     */
    void setNx(const Real Nx) { M_Nx = Nx; }

    //! Set the number of cells along y-axis
    /*!
     * @param the number of cells along y-axis
     */
    void setNy(const Real Ny) { M_Ny = Ny; }

    //! Set the number of cells along z-axis
    /*!
     * @param the number of cells along z-axis
     */
    void setNz(const Real Nz) { M_Nz = Nz; }

    //! Set the shift domain along x-axis
    /*!
     * @param the shift domain along x-axis
     */
    void setSx(const Real Sx) { M_Sx = Sx; }

    //! Set the shift domain along y-axis
    /*!
     * @param the shift domain along y-axis
     */
    void setSy(const Real Sy) { M_Sy = Sy; }

    //! Set the shift domain along z-axis
    /*!
     * @param the shift domain along z-axis
     */
    void setSz(const Real Sz) { M_Sz = Sz; }

    //! Set the numerical method
    /*!
     * @param type the type of the numerical method
     */
    void setNumericalMethodType(const NumericalMethodType type) { M_numet = type; }

    //! Set if the stiffness matrix in the mimetic method is lumped
    /*!
     * @param lumped true if the stiffness matrix is lumped
     */
    void setLumpedMimetic(const bool lumped) { M_lumpedMim = lumped; }

    //! Set the problem type
    /*!
     * @param type the type of the problem
     */
    void setProblemType(const ProblemType type) { M_problemType = type; }

    //! Enable or disable the fractures
    /*!
     * @param fracture true to enable the fractures
     */
    void fractureOn(bool fracture) { M_fracturesOn = fracture; }

    //! Set where the source/sink term is applied
    /*!
     * @param ssOn where the source/sink term is applied
     */
    void setSourceSinkOn(const SourceSinkOn ssOn) { M_ssOn = ssOn; }

    //! Enable or disable the pressure value inside the fractures
    /*!
     * @param fracPress if true fix a pressure value inside the fractures
     */
    void pressuresInFractures(bool fracPress) { M_setFracturesPressure = fracPress; }

    //! Set the pressure inside the fracture
    /*!
     * @param press the pressure inside the fracture
     */
    void setPressuresInFractures(const Real press) { M_fracturesPressure = press; }

    //! Enable or disable the the multiple sub-regions method
    /*!
     * @param msr if true the multiple sub-regions method is activated
     */
    void MSROn(bool msr) { M_MSR = msr; }

    //! Set the number of sub-regions
    /*!
     * @param nSubReg the number of sub-regions
     */
    void setNbSubRegions(const UInt nSubReg) { M_nbSubRegions = nSubReg; }

    //! Set the number of time step to check for the steady state
    /*!
     * @param nTimeStep the number of time step to check for the steady state
     */
    void setNbTimeStepSteadyState(const UInt nTimeStep ) { M_nbTimeStepSteadyState = nTimeStep; }

    //! Set the tolerance to consider a time step at the steady state
    /*!
     * @param tol the tolerance to consider a time step at the steady state
     */
    void tolSteadyState(const Real tol) { M_tolSteadyState = tol; }

    //! Set the type of the permeability
    /*!
     * @param type the type of the permeability
     */
    void setPermeabilityType(const std::string type) { M_permeabilityType = type; }

    //! Set the permeability in the porous medium
    /*!
     * @param perm the permeability in the porous medium
     */
    void setMatrixPermeability(const Real perm) { M_permMatrix = perm; }

    //! Set the porosity in the porous medium
    /*!
     * @param poro the porosity in the porous medium
     */
    void setMatrixPorosity(const Real poro) { M_poroMatrix = poro; }

    //! Set the permeability in the fracture
    /*!
     * @param perm the permeability in the fracture
     */
    void setFracturePermeability(const Real perm) { M_permFrac = perm; }

    //! Set the porosity in the fracture
    /*!
     * @param poro the porosity in the fracture
     */
    void setFracturePorosity(const Real poro) { M_poroFrac = poro; }

    //! Set the aperture in the fracture
    /*!
     * @param aper the aperture in the fracture
     */
    void setFractureAperture(const Real aper) { M_aperFrac = aper; }

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

    //! Set solver type
    /*!
     * @param the solver type
     */
    void getSolverType(const std::string solver) { M_solverType = solver; }

    //! Set maximum iterations of the iterative solver
    /*!
     * @param maxIt maximum iterations of the iterative solver
     */
     void setIterativeSolverMaxIter(const UInt maxIt) { M_maxIt = maxIt; }

    //! Set tolerance of the iterative solver
    /*!
     * @param tol tolerance of the iterative solver
     */
    void setIterativeSolverTolerance(const Real tol) { M_tol = tol; }

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
    //! Domain dimension along x-axis
    Real M_Lx;
    //! Domain dimension along y-axis
    Real M_Ly;
    //! Domain dimension along z-axis
    Real M_Lz;
    //! Number of cells along x-axis
    UInt M_Nx;
    //! Number of cells along y-axis
    UInt M_Ny;
    //! Number of cells along z-axis
    UInt M_Nz;
    //! Shift domain along x-axis
    Real M_Sx;
    //! Shift Domain along y-axis
    Real M_Sy;
    //! Shift Domain along z-axis
    Real M_Sz;
    //! Type of the numerical method
    NumericalMethodType M_numet;
    //! If true the stiffness matrix for the mimetic method is lumped
    bool M_lumpedMim;
    //! Type of the problem
    ProblemType M_problemType;
    //! Enable or disable fractures
    bool M_fracturesOn;
    //! Where the source/sink term is applied
    SourceSinkOn M_ssOn;
    //! If true fix the pressure inside the fractures
    bool M_setFracturesPressure;
    //! Pressure value inside the fractures
    Real M_fracturesPressure;
    //! Multiple sub-regions on?
    bool M_MSR;
    //! Number of sub-regions
    UInt M_nbSubRegions;
    //! Number of time step to check for the steady state
    UInt M_nbTimeStepSteadyState;
    //! Tolerance to consider a time step at the steady state
    Real M_tolSteadyState;
    //! Type of the permeability
    std::string M_permeabilityType;
    //! Permeability in the porous medium
    Real M_permMatrix;
    //! Porosity in the porous medium
    Real M_poroMatrix;
    //! Permeability in the fractures
    Real M_permFrac;
    //! Porosity in the fractures
    Real M_poroFrac;
    //! Aperture in the fractures
    Real M_aperFrac;
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

    //! Solver used
    std::string M_solverType;
    //! Maximum iterations of the iterative solver
    UInt M_maxIt;
    //! Tolerance of the iterative solver
    Real M_tol;

    //! Angle used to rotate along z-axis the domain. It is used only to compute the normal for detecting BC!
    Real M_theta;
    //! Verbose
    bool M_verbose;
};

//! Class that implements a parser for a generic enumerator
/*!
 * @class EnumParser
 * This class allows to convert a string to an enumerator value
 * The template parameter is the enumerator type
 */
template <class T>
class EnumParser
{
public:

    //! Empty constructor
    EnumParser();

    //! Parse method
    /*!
     * @param str string of the enumerator value that you want to convert
     * @return the enumerator value associated with the string
     */
    T parse(const std::string & str) const throw()
    {
        typename std::map<std::string, T>::const_iterator it = M_enumMap.find(str);
        if (it == M_enumMap.end())
        {
            throw std::runtime_error("Error: parsing an enum that does not exist.");
        }
        return it->second;
    }

    //! Destructor
    ~EnumParser() = default;

private:

    //! Map that links a string to an enumarator type
    std::map<std::string, T> M_enumMap;
};

//! Specialized class that implements a parser for the Data::NumericalMethodType enum
template<>
EnumParser<Data::NumericalMethodType>::EnumParser();

//! Specialized class that implements a parser for the Data::ProblemType enum
template<>
EnumParser<Data::ProblemType>::EnumParser();

//! Specialized class that implements a parser for the Data::MeshFormatType enum
template<>
EnumParser<Data::MeshFormatType>::EnumParser();

//! Specialized class that implements a parser for the Data::SourceSinkOn enum
template<>
EnumParser<Data::SourceSinkOn>::EnumParser();

//! Typedef for std::unique_ptr<Data>
/*!
 * @typedef DataPtr_Type
 * This type definition permits to handle a std::unique_ptr<Data> as a DataPtr_Type.
 */
typedef std::unique_ptr<Data> DataPtr_Type;

} // namespace FVCode3D
#endif /* DATA_HPP_ */
