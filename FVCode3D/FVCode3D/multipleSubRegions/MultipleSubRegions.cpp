/*!
 *  @file MultipleSubRegions.cpp
 *  @brief Classes that implements the multiple sub-regions (definitions).
 */

#include "multipleSubRegions/MultipleSubRegions.hpp"

namespace FVCode3D
{

UInt linearRegionSelection ( const Real& _cellSolution, const UInt& _nbRegions,
                             const Real& _minSolution, const Real& _maxSolution )
{
    const Real tol = 10. * std::numeric_limits<Real>::epsilon();
    const Real newMinSolution = _minSolution * ( 1. + tol * ( (_minSolution > 0.) ? -1. :  1. ) );
    const Real newMaxSolution = _maxSolution * ( 1. + tol * ( (_maxSolution > 0.) ?  1. : -1. ) );
    return ceil ( _nbRegions * std::fabs( _cellSolution - newMinSolution  ) / ( std::fabs( newMaxSolution - newMinSolution ) ) );
} // linearRegionSelection

} // namespace FVCode3D
