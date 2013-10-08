#ifndef HH_TRILINEARSHAPE_HH
#define HH_TRILINEARSHAPE_HH
#include <stdexcept>

namespace Geometry
{

  /*! Shape functions for linear elements.
    
    End nodes are numbered 0 and 1 and passed as template argument.
  */

  template <unsigned short I>
  inline double linearShape (double const& x);

  //! Full specialization.
  template<>
  inline double linearShape<0> (double const& x)
  {
    return (1.0 - x);
  }

  //! Full specialization.
  template<>
  inline double linearShape<1> (double const& x)
  {
    return x;
  }

  /*!
    \f$ \partial\psi_I/\partial x\f$.
   */
  template<unsigned short I>
  inline const double  linearShapeDer (double const& x);

  template<>
  inline const double linearShapeDer<0> (double const& x)
  {
    return -1.0;
  }

  template<>
  inline const double linearShapeDer<1> (double const& x)
  {
    return 1.0;
  }

  /*! Reference shape functions for a trilinear finite element.

    We use the following node numbering

\verbatim
         6 ________________.7
          /.              /!
         / .             / !
       4!---------------!5 !
        !  .            !  !
        !  .            !  !
        !  .2...........! .!3
        !.              ! /
      0 _________________/1

\endverbatim

      and its binary representation K-> K mod(2),K/2 mod 2, K/4 mod 2.

      Node numbering (in its decimal representation) is passed as
      template argument. Its binary representation is used to select
      the appropriate linear shape functions (since a trilinear
      function is the product of three linear shape functions).
  */
  template<unsigned short I>
  inline double trilinearShape (double const& x, double const& y, double const& z)
  {
    return linearShape < I % 2 > (x) * linearShape < (I / 2) % 2 > (y) * linearShape < (I / 4) % 2 > (z);
  }

  /*!\f$ \partial\phi_I/\partial x_J\f$.
   */
  template<unsigned short I, unsigned short J>
  inline double trilinearShapeDer (double const& x, double const& y, double const& z)
  {
    switch (J)
    {
      case 0:
        return linearShapeDer < I % 2 > (x) * linearShape < (I / 2) % 2 > (y) *
               linearShape < (I / 4) % 2 > (z);
        break;
      case 1:
        return linearShape < I % 2 > (x) * linearShapeDer < (I / 2) % 2 > (y) *
               linearShape < (I / 4) % 2 > (z);
        break;
      case 2:
        return linearShape < I % 2 > (x) * linearShape < (I / 2) % 2 > (y) *
               linearShapeDer < (I / 4) % 2 > (z);
        break;
      default:
        throw std::range_error ("Coordinate index is either 0, 1 or 2");
    }
  }

}// end namespace Geometry
#endif
