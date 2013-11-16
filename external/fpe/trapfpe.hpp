#ifndef __TRAPFPE_HPP__
#define __TRAPFPE_HPP__
#ifdef FPE_ABORT
#include <cfenv>
/*! \file trapfpe.hpp An example of use of floating point environment

  We have a cpp macro
  FPE_ABORT Activates abort on floating point exception
  
  If only FPE_ON is activated just the function test_fpe_exception
  is compiled.
  
  You should include this file anc compile with -DFPE_ABORT if you want to 
  abort on a floating point exception.
 */
// The C++11 standard requires this pragma to be activated
// Only some compilers have this feature however. Ignore the warning
#pragma STDC FENV_ACCESS ON
#warning "ABORTING ON FPE IS ACTIVATED" 
//! If FPE_ABORT is set a floatin point exception causes the progam to terminate
/* C style implementation
   static void __attribute__ ((constructor))
   trapfpe (){
   // Enable some exceptions.  At startup all exceptions are masked.
feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}
*/
//! IMPLEMENTATION AS STATIC METHOD (More C++ style)

namespace{
  //! Structure that wraps the function enabling the prapping of fpe exceptions.
  struct FpeTrap{
      /*! \brief Enable some exceptions.  

	\detail We consider only the following FPEs:
	Invalid operation (i.e. log(-1)), division by zero and overflow.
      */
    static void __attribute__ ((constructor))
    trapfpe (){
      feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
    }
  };
}
#endif
#endif
