/*
 *  @file data.hpp
 *  @brief This class handles the data for the solver.
 */

#ifndef DATA_HPP_
#define DATA_HPP_

#include <FVCode3D/core/TypeDefinition.hpp>
#include <FVCode3D/utility/StringManipolator.hpp>
#include <FVCode3D/utility/readPermeabilities.hpp>
#include <FVCode3D/utility/readOtherData.hpp>

namespace FVCode3D
{
  
  //! Class that handles the information needed for the solver
  /*!
   * @class Data
   * This class takes care of the information needed for the solver, as boundary conditions and the source/sink term.
   * It also handles the name of the input/output directory and the name of the input/output file,
   * the shared library name and miscellaneous options.
   * @note This class is a monster class overly complicated. Its role is to have a centralized place for all problem parameters.
   * The idea is fine yet it had been better to have a simple aggregate (no need of setters and getters) maybe implemented as
   * a singleton namespace variable (no need of using pointers in this case). It would make life much easier.
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
    
    //! Define the solver policy
    /*!
     * @enum SolverPolicy
     * This enumerator allows to select the solver policy: direct solver or iterative solver
     */
    enum SolverPolicy
      {
        Direct     = 0,
        Iterative  = 1
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
     * This enumerator allows to select the format type of the mesh:
     * TPFA(.grid), forSolver(.fvg), Medit(.mesh) or TetGen(.node, .face, .ele)
     */
    enum MeshFormatType
      {
        TPFA            = 0,
        forSolver       = 1,
        Medit           = 2,
        TetGen          = 3,
        OpenFOAM        = 4
      };
    
    //! Define where the noise is applied
    /*!
     * @enum NoiseOn
     * This enumerator allows to select where the noise is applied:
     * only on matrix nodes, only on fracture nodes, on all nodes.
     */
    enum class NoiseOn
    {
      Matrix          = 0,
      Fractures       = 1,
      All             = 2
      };
    enum BcStrategy
    {
      Strong = 0,
      Nitsche= 1
    };
    //! Define where to apply the source/sink term
    /*!
     * @enum SourceSinkOn
     * This enumerator allows to select where to apply the source/sink term: only on matrix, only on fractures, both or none.
     */
    enum class SourceSinkOn
    {
      Matrix          = 0,
      Fractures       = 1,
      Both            = 2,
      None            = 3
    };

    //! Defines permeability types
    /*!
     * @enum PermeabilityType
     * Allows to select the type of permeability
     * \note ArbirtraryTensor is kept for compatiblity, but it is not used in practice so far.
     */
     enum class PermeabilityType
     {
   	  Scalar=0,
   	  Diagonal=1,
   	  SymTensor=2,
   	  ArbitraryTensor=3
     };

    //! Defines the inner solver used for the Schur complement in the preconditioner
     enum class SchurInnerSolverType
     {
       Chowlesky=0,
       Umfpack=1,
       CG=2
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

    //! Load data from filename
    void load(const std::string dataFileName);
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

    //! Get the rotation domain z-axis
    /*!
     * @return the rotation domain along z-axis
     */
    Real getRz() const { return M_Rz; }

    //! Test if the noise is added to the points
    /*!
     * @return true if the noise is added to the points
     */
    bool noiseOn() const { return M_noise; }

    //! Get where the noise is applied
    /*!
     * @return where the noise is applied
     */
    NoiseOn getNoiseOn() const { return M_noiseOn; }

    //! Get the mean of the normal distribution
    /*!
     * @return the mean of the normal distribution
     */
    Real getMeanNormalDistribution() const { return M_meanNDist; }

    //! Get the standard deviation of the normal distribution
    /*!
     * @return the standard deviation of the normal distribution
     */
    Real getStDevNormalDistribution() const { return M_stDevNDist; }

    //! Get the numerical method used to solve the problem
    /*!
     * @return the numerical method used to solve the problem
     */
    NumericalMethodType getNumericalMethodType() const { return M_numet; }

    //! Get if the stiffness matrix in the mimetic method is lumped
    /*!
     * @return lumped true if the stiffness matrix is lumped
     * This method is never called. We need to fix it!
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

    //! Get matrix (bulk) permeability associated to a zone.
    /*!
     * @par zoneNumber The number of the zone. 0 default
     * @return A tuple<int,vector<Real>> containing the type and the values
     * The type is according to PermeabilityType, the values are accordingly
     *
     * @throw runtime_error If zone is not found
     */
    Utility::Permeability getMatrixPermeability(int zoneNumber) const;
    //! Get the permeability in the porous medium
     /*!
      * @return the permeability in the porous medium
      */
    Real getMatrixPermeability() const { return M_permMatrix; }
    //! Get the permeability in the porous medium as a PermeabilityDataList
    /*!
     * @return the permeability in the porous medium
     */
    Utility::BulkPermeabilityData getMatrixPermeabilityData() const
    {
      return this->M_bulkPermeabilityData;
    }
    //! Get the porosity in the porous medium
    /*!
     * @return the porosity in the porous medium
     */
    Real getMatrixPorosity() const { return M_poroMatrix; }
    //! Get the porosity in the porous medium
    /*!
     * @return the porosity in the porous medium as a bulkDataList
     */
    Utility::bulkDataList getMatrixOtherData() const { return this->M_bulkData; }
    //! Get matrix (bulk) porosity associated to a zone.
    /*!
     * @par zoneNumber The number of the zone. 0 default
     * @return A Real containing the corresponding value     *
     * @throw runtime_error If zone is not found
     */
    Real getMatrixPorosity(int zoneNumber) const;

    //! Get fracture permeability associated to a zone.
    /*!
     * @par zoneNumber The number of the zone. 0 default
     * @return A tuple<int,vector<Real>> containing the type and the values
     * The type is according to PermeabilityType, the values are accordingly
     *
     * @throw runtime_error If zone is not found
     */
    Utility::Permeability getFracturePermeability(int zoneNumber) const;
    //! Get the permeabilities in the various zones of the fracture as a PeremabilityDataList item
    /*!
     * @return the permeability in the fracture
     */
    Utility::FracturePermeabilityData getFracturePermeabilityData() const
    {
      return this->M_fracturePermeabilityData;
    }
    //! Get the permeability in the fracture as scalar
    /*!
     * @return the permeability in the fracture
     */

    Real getFracturePermeability() const { return M_permFrac; }

    //! Get fracture porosity associated to a zone.
    /*!
     * @par zoneNumber The number of the zone. 0 default
     * @return A Real containing the corresponding value     *
     * @throw runtime_error If zone is not found
     */
    Real getFracturePorosity(int zoneNumber) const;

    //! Get fracture porosity and aperture for each zone in the form of a FractureData List
      /*!
     * @return A FractureDataList object containing data on porosity and aperture
     * fo each fracture zone
     */
    Utility::fractureDataList getFractureOtherData() const
    {
      return this->M_fractureData;
    }

    //! Get the porosity in the fracture
    /*!
     * @return the porosity in the fracture
     */
    Real getFracturePorosity() const { return M_poroFrac; }

    //! Get fracture aperture associated to a zone.
    /*!
     * @par zoneNumber The number of the zone. 0 default
     * @return A Real containing the corresponding value     *
     * @throw runtime_error If zone is not found
     */
    Real getFractureAperture(int zoneNumber) const;

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

    //! Get solver policy
    /*!
     * @return the solver policy
     */
    SolverPolicy getSolverPolicy() const { return M_SolverPolicy; }

    //! Get solver type
    /*!
     * @return the solver type
     */
    const std::string getSolverType() const { return M_solverType; }
    
    //! Get preconditioner type
    /*!
     * @return the preconditioner type
     */
    const std::string getpreconType() const { return M_precon; }

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

    //! Get restart (only for GMRES and FGMRES)
    //! @return the restart for GMRES(m) and FGMRES(m)
    //! @todo This method is never used, I need to fix it.
     UInt getRestart(){return M_restart;}

    //! Get theta angle
    /*!
     * @return the angle theta used to identify the BC ids
     */
    Real getTheta() const { return M_theta; }

    //! Do we want to dump matrices in matrixmarket format?
    /*
     * @return true. System matrices are dumped in matrix market format with standard names
     */
    bool getDumpMatrices()const {return M_dumpMatrix;}
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

    //! Set the rotation domain along z-axis
    /*!
     * @param the rotation domain along z-axis
     */
    void setRz(const Real Rz) { M_Rz = Rz; }

    //! Enable or disable the noise to the points
    /*!
     * @param noise true to enable the noise to the points
     */
    void noiseOn(const bool noise) { M_noise = noise; }

    //! Set where the noise is applied
    /*!
     * @param noiseOn where the noise is applied
     */
    void getNoiseOn(const NoiseOn noiseOn) { M_noiseOn = noiseOn; }

    //! Set the mean of the normal distribution
    /*!
     * @param mean the mean of the normal distribution
     */
    void setMeanNormalDistribution(const Real mean) { M_meanNDist = mean; }

    //! set the standard deviation of the normal distribution
    /*!
     * @param stDev the standard deviation of the normal distribution
     */
    void setStDevNormalDistribution(const Real stDev) { M_stDevNDist = stDev; }

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
    void fractureOn(const bool fracture) { M_fracturesOn = fracture; }

    //! Set where the source/sink term is applied
    /*!
     * @param ssOn where the source/sink term is applied
     */
    void setSourceSinkOn(const SourceSinkOn ssOn) { M_ssOn = ssOn; }

    //! Enable or disable the pressure value inside the fractures
    /*!
     * @param fracPress if true fix a pressure value inside the fractures
     */
    void pressuresInFractures(const bool fracPress) { M_setFracturesPressure = fracPress; }

    //! Set the pressure inside the fracture
    /*!
     * @param press the pressure inside the fracture
     */
    void setPressuresInFractures(const Real press) { M_fracturesPressure = press; }

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

    //! Set solver policy
    /*!
     * @param the solver policy
     */
    void setSolverPolicy(const SolverPolicy solver) { M_SolverPolicy = solver; }

    //! Set preconditioner type
    /*!
     * @param the preconditioner type
     */
    void setpreconType(const std::string precon) { M_precon =precon; }
    //! Set inner solver type for the approximate Schur complement
    void setSchurInnerSolverType(SchurInnerSolverType s){this->M_SchurInnerSolverType = s;}
    //! Set tolerance for inner solver for the Schur complement, if iterative
    void setSchurInnerSolverTolerance(Real s){this->M_SchurInnerSolverTolerance = s;}
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

    //! Set restart (only for GMRES and FGMRES)
    void setRestart(UInt r){this->M_restart=r;}
    //! Set theta angle
    /*!
     * @param theta the angle theta used to identify the BC ids
     */
    void setTheta(const Real theta) { M_theta = theta; }
    //! Sets if we want to dump matrices in matrixmarket format
    void setDumpMatrices(bool dump){this->M_dumpMatrix=dump;}

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

public:

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
    //! Rotation Domain along z-axis
    Real M_Rz;
    //! If true noise is added to the points
    bool M_noise;
    //! Where the noise is applied
    NoiseOn M_noiseOn;
    //! Mean of the normal distribution
    Real M_meanNDist;
    //! Standard deviation of the normal distribution
    Real M_stDevNDist;
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

    //! Solver policy
    SolverPolicy M_SolverPolicy;
    //! Solver used
    std::string M_solverType;
    //! Maximum iterations of the iterative solver
    UInt M_maxIt;
    //! Tolerance of the iterative solver
    Real M_tol;
    //! restart (for GMRES and FGMRES)
    UInt M_restart=100u;
    //! The preconditioner
    std::string M_precon;
    //! Lumping or Diagonal?
    // bool M_lumped; // now M_lumpedMim
    //! InnerIterativeSolver
    SchurInnerSolverType M_SchurInnerSolverType={SchurInnerSolverType::Chowlesky};
    //! Tolerance for the Inner Iterative Solver
    Real M_SchurInnerSolverTolerance={1.e-7};
    //! Angle used to rotate along z-axis the domain. It is used only to compute the normal for detecting BC!
    Real M_theta;
    //! Verbose
    bool M_verbose;
    //      Version 2 by Luca F.
    //! Contains data on bulk permeability for each zone
    Utility::BulkPermeabilityData M_bulkPermeabilityData;
    //! Contains data on bulk permeability for each zone
    Utility::FracturePermeabilityData M_fracturePermeabilityData;
    //! Contains porosity and aperture for each zone
    Utility::bulkDataList M_bulkData;
    //! Contains porosity for each zone
    Utility::fractureDataList M_fractureData;
    //! Do we want to dump matrices?
    /*!
     * If true all matrices are dumped in matrixmarket format
     */
    bool M_dumpMatrix{false};
    //! Strategy for bc imposition
    BcStrategy M_bcStrategy=Nitsche;
    //! Nitsche penalization parameter
    double M_Nitsche=100.;
    //
private:
    //! Utility to read fracture porosity and aperture data
    std::tuple<double,double> fracturePorosityAndAperture(int const & zoneNumber) const;
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
    
    //! Default constructor
    EnumParser();
    
    //! Parse method
    /*!
     * @param str string of the enumerator value that you want to convert
     * @return the enumerator value associated with the string
     */
    T parse(const std::string & str) const
    {
      typename std::map<std::string, T>::const_iterator it = M_enumMap.find( toUpper( str ) );
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

//! Specialized class that implements a parser for the Data::SolverPolicy enum
template<>
EnumParser<Data::SolverPolicy>::EnumParser();

//! Specialized class that implements a parser for the Data::ProblemType enum
template<>
EnumParser<Data::ProblemType>::EnumParser();

//! Specialized class that implements a parser for the Data::MeshFormatType enum
template<>
EnumParser<Data::MeshFormatType>::EnumParser();

//! Specialized class that implements a parser for the Data::NoiseOn enum
template<>
EnumParser<Data::NoiseOn>::EnumParser();

//! Specialized class that implements a parser for the Data::SourceSinkOn enum
template<>
EnumParser<Data::SourceSinkOn>::EnumParser();

//! Specialized class that implements a parser for the Data::SourceSinkOn enum
template<>
EnumParser<Data::PermeabilityType>::EnumParser();

//! Typedef for std::unique_ptr<Data>
/*!
 * @typedef DataPtr_Type
 * This type definition permits to handle a std::unique_ptr<Data> as a DataPtr_Type.
 */
typedef std::unique_ptr<Data> DataPtr_Type;

//! This is a new treatment for Data object. We access them through a namespace global object

extern DataPtr_Type dataPtr;

} // namespace FVCode3D


#endif /* DATA_HPP_ */
