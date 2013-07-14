#ifndef HH__TOLERANCES_HH
#define HH__TOLERANCES_HH
#include "TypeDefinition.hpp"
//! A namespace holding the tolerances used in the code
/*!
  @author Luca Formaggia

  We want to keep tolerances well enucleated so that it is easier to modify
  and check them. This namespace serves this purpose.
  To make cleaar that they are constants we write them all uppercase.

  @note If we were using C++11 we could declare the tolerances as constexpr. 
  
 */
namespace EDFM_Tolerances
{
  //! A base tolerance.
  Real const BASE_TOLERANCE=1.e-8;
  //! Used in Newton algorithms.
  Real const NEWTON_TOLERANCE=1.e-6;
  //! For alignement of points
  Real const ALIGNMENT_TOLERANCE=1.e-5;
}// end namespace EDFM_Tolerances
#endif
