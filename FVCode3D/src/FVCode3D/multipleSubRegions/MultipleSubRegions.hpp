 /*!
 * @file MultipleSubRegions.hpp
 * @brief Classes that implements the multiple sub-regions.
 */

#ifndef MULTIPLESUBREGIONS_HPP_
#define MULTIPLESUBREGIONS_HPP_

#include <functional>
#include <memory>

#include <FVCode3D/eigenPatch/RangeSupport.hpp>
#include <FVCode3D/core/TypeDefinition.hpp>
#include <FVCode3D/core/Data.hpp>

namespace FVCode3D
{

//! Given a cell solution, this function linearly divides the solution range and returns in which sub-region the cell belongs
/*!
 * @param _cellSolution solution in a given cell
 * @param _nbRegions number of sub regions
 * @param _minSolution the minimum value assumed by the solution
 * @param _maxSolution the maximum value assumed by the solution
 * @return the sub-region id that cell belongs
 */
UInt linearRegionSelection ( const Real& _cellSolution, const UInt& _nbRegions,
                             const Real& _minSolution = 0, const Real& _maxSolution = 0 );

//! Class that implements the multiple sub regions method
/*!
 * @class MSR
 * This class implements the multiple sub regions method.
 * Given the solution and the number of sub regions, it computes for each cell in which sub region the cell belongs.
*/
template <typename ProblemType>
class MSR
{
public:

    typedef ProblemType problem_Type;
    typedef MSR<problem_Type> class_Type;
    typedef std::function<UInt(const Real&, const UInt&, const Real&, const Real& )> regionFct_Type;

    //! @name Constructors and destructor
    //@{

    //! No empty constructor
    MSR() = delete;

    //! Constructor
    /*!
     * @param problem a pointer of a Problem class
     * @param data a constant reference to a Data class
     */
    MSR(ProblemType* problem, const DataPtr_Type& data);

    //! Default destructor
    ~MSR() = default;

    //@}

    //! @name Get methods
    //@{

    //! Get the vector with the id of the sub-region for each cell
    /*!
     * @return the vector with the id of the sub-region for each cell
     */
    const UIntVector& getColorsVector() const { return M_multipleSubRegions; }

    //! Get the transmissibilities vector
    /*!
     * @return the transmissibilities vector: fracture-matrix and matrix-matrix transmissimilities
     */
    const std::vector<Real>& getTransmissibility() const { return M_transmissibilities; }

    //@}

    //! @name Methods
    //@{

    //! Set dofs and initialize the sub-regions vector
    void setup();

    //! Compute for each cell in which region it belongs
    /*!
     * @param _regionFct function used to divide the range solution
     */
    void createRegions(const regionFct_Type& _regionFct = linearRegionSelection);

    //! Check if the solution is at the pseudo steady state
    /*!
     * @param _oldSolution solution at the previous time step
     */
    bool isPseudoSteadyState(const Vector& _oldSolution);

    //! Compute the transmissibility between each sub region
    /*!
     * It also compute the transmissibility at the domain boundary
     */
    void computeTransmissibility();

    //@}

private:

    //! Compute the variation of the pressure
    inline void computeDeltaPressure();

    //! Compute the average pressure in each sub region
    /*!
     * @param _averagePressure vector that will contain the average pressures
     * The last value contains the boundary pressure
     */
    void computeAveragePressure(std::vector<Real>& _averagePressure) const;

    //! Compute the flux between each sub region
    /*!
     * @param flux vector that will contain the flux
     * It also compute the transmissibility at the domain boundary
     */
    void computeFlux(std::vector<Real>& flux) const;

    //! Darcy problem
    problem_Type* M_problem;
    //! Data for the problem
    const DataPtr_Type& M_data;

    //! Problem total dofs
    UInt M_dofTotal;
    //! Fracture dofs
    UInt M_dofFractures;
    //! Porous matrix dofs
    UInt M_dofPorousMatrix;
    //! Vector of colors of the sub-regions
    UIntVector M_multipleSubRegions;
    //! Maximum of the solution
    Real M_maxSolution;
    //! Minimum of the solution
    Real M_minSolution;
    //! If the state is a pseudo-steady state
    bool M_isPseudoSteadyState;
    //! Number of the time step that the state is at the pseudo-steady state
    UInt M_localNbTimeStepSteadyState;
    //! c := \frac{ \partial p }{ \partial t}
    Real M_constantValue;
    //! Transmissibilities vector
    std::vector<Real> M_transmissibilities;
}; // class MSR

template <typename ProblemType>
MSR<ProblemType>::MSR(ProblemType* problem, const DataPtr_Type& data ):
    M_problem( problem ),
    M_data( data ),
    M_dofTotal( 0 ),
    M_dofFractures( 0 ),
    M_dofPorousMatrix( 0 ),
    M_maxSolution( std::numeric_limits<Real>::lowest() ),
    M_minSolution( std::numeric_limits<Real>::max() ),
    M_isPseudoSteadyState(false),
    M_localNbTimeStepSteadyState(0),
    M_constantValue( std::numeric_limits<Real>::lowest() )
{} // MSR::MSR

template <typename ProblemType>
void MSR<ProblemType>::setup()
{
    M_dofTotal = M_problem->getA().outerSize();
    M_dofFractures = M_problem->getMesh().getFractureFacetsIdsVector().size();
    M_dofPorousMatrix = M_dofTotal - M_dofFractures;
    M_multipleSubRegions = UIntVector::Zero( M_dofTotal );
} // MSR::setup

template <typename ProblemType>
void MSR<ProblemType>::createRegions( const regionFct_Type& _regionFct )
{
    computeDeltaPressure();

    using namespace std::placeholders;
    const UInt nbSubRegions = M_data->nbSubRegions();
    const auto fc_region = std::bind( _regionFct, _1, nbSubRegions, M_minSolution, M_maxSolution );

    size_t index = 0;
    const auto & solution = M_problem->getSolver().getSolution();
    auto apply_sub_region = [&] ( UInt& value ){ value = fc_region( solution[ index ] ); ++index; };

    std::for_each( std::begin( M_multipleSubRegions ), std::end( M_multipleSubRegions ), apply_sub_region );
} // MSR::createRegions

template <typename ProblemType>
inline void MSR<ProblemType>::computeDeltaPressure()
{
    const auto & solution = M_problem->getSolver().getSolution();

    M_maxSolution = *std::max_element( std::begin( solution ), std::end( solution ) );
    M_minSolution = *std::min_element( std::begin( solution ), std::end( solution ) );
} // MSR::computeDeltaPressure

template <typename ProblemType>
bool MSR<ProblemType>::isPseudoSteadyState(const Vector& _oldSolution )
{
    const Vector& solution = M_problem->getSolver().getSolution();
    const Real currentConstantValue = ( solution - _oldSolution ).norm() / M_data->getTimeStep();
    M_isPseudoSteadyState = ( currentConstantValue - M_constantValue ) < M_data->tolSteadyState();
    M_constantValue = currentConstantValue;
    M_localNbTimeStepSteadyState++;

    std::cout << " currentConstantValue " << currentConstantValue << std::endl;
    std::cout << " diff " << currentConstantValue - M_constantValue << std::endl;
    std::cout << " M_isPseudoSteadyState " << M_isPseudoSteadyState << std::endl;

    return ( M_isPseudoSteadyState && M_localNbTimeStepSteadyState >= M_data->nbTimeStepSteadyState() );
} // MSR::isPseudoSteadyState

template <typename ProblemType>
void MSR<ProblemType>::computeTransmissibility()
{
    const UInt nbSubRegions = M_data->nbSubRegions();
    std::vector<Real> averagePressure( nbSubRegions + 1, 0. ), flux( nbSubRegions, 0. );
    computeAveragePressure( averagePressure );
    computeFlux( flux );

    M_transmissibilities.resize( nbSubRegions );
    for ( size_t i = 0; i < M_transmissibilities.size(); ++i )
    {
        const Real deltaPressure = std::fabs ( averagePressure[ i + 1 ] - averagePressure[ i ] );
        M_transmissibilities[ i ] = flux[ i ] / deltaPressure;
    } // for
} // MSR::computeTransmissibility

template <typename ProblemType>
void MSR<ProblemType>::computeAveragePressure( std::vector<Real> & _averagePressure ) const
{
    const Vector& solution = M_problem->getSolver().getSolution();
    const auto& properties = M_problem->getMesh().getPropertiesMap();
    const auto& cellsVector = M_problem->getMesh().getCellsVector();
    const auto& cellsFault = M_problem->getMesh().getFractureFacetsIdsVector();
    std::vector<Real> measures( M_multipleSubRegions.size(), 0. );
    std::vector<Real> area( _averagePressure.size(), 0. );

    for(auto& cell_it : cellsVector )
    {
        measures[ cell_it.getId() ] = cell_it.getVolume();
    } // for

    for(auto& facet_it : cellsFault )
    {
        const Real localMeasure = properties.getProperties( facet_it.getZoneCode() ).M_aperture
                                  * facet_it.getFacet().area();
        measures[ facet_it.getIdAsCell() ] = localMeasure;
        const auto& neighborCells = facet_it.getSeparatedCellsIds();
        for(auto& cell_it : neighborCells )
        {
            measures[ cell_it ] -= localMeasure / 2.;
        } // for
    } // for

    for(size_t i = 0; i < M_multipleSubRegions.size(); ++i )
    {
        const UInt regionId = M_multipleSubRegions[ i ];
        const Real cellMeasure = measures[ i ];
        _averagePressure[ regionId ] += solution[ i ] * cellMeasure;
        area[ regionId ] += cellMeasure;
    } // for

    for(size_t i = 0; i < _averagePressure.size(); ++i )
    {
        _averagePressure[ i ] /= area[ i ];
    } // for
} // MSR::computeAveragePressure

template <typename ProblemType>
void MSR<ProblemType>::computeFlux( std::vector<Real> & _flux ) const
{
    const auto& solution = M_problem->getSolver().getSolution();
    const auto& mesh = M_problem->getMesh();
    for(auto& facet_it : mesh.getInternalFacetsIdsVector() )
    {
        const UInt neighborCell0 = facet_it.getSeparatedCellsIds()[0];
        const UInt neighborCell1 = facet_it.getSeparatedCellsIds()[1];
        const UInt colorCell0 = M_multipleSubRegions[ neighborCell0 ];
        const UInt colorCell1 = M_multipleSubRegions[ neighborCell1 ];
        if ( colorCell0 != colorCell1 )
        {
            assert( std::abs( static_cast<Int>( colorCell0 - colorCell1 ) ) == 1 );
            const Real deltaSolution = std::fabs ( solution[ neighborCell0 ] - solution[ neighborCell1 ] );
            const Real transmissibility = M_problem->getA().coeffRef( neighborCell0, neighborCell1 );
            _flux[ std::min( colorCell0, colorCell1 ) ] += transmissibility * deltaSolution;
        } // if
    } // for

    for(auto& facet_it : mesh.getFractureFacetsIdsVector() )
    {
        const auto& neighborCells = facet_it.getSeparatedCellsIds();
        const UInt facetIdasCell = facet_it.getIdAsCell();
        for(auto& cell_it : neighborCells )
        {
            const Real transmissibility = M_problem->getA().coeffRef( cell_it, facetIdasCell );
            const Real deltaSolution = std::fabs ( solution[ cell_it ] - solution[ facetIdasCell ] );
            _flux[ 0 ] += transmissibility * deltaSolution;
        } // for
    } // for
} // MSR::computeFlux

} // namespace FVCode3D

#endif /* MULTIPLESUBREGIONS_HPP_ */
