/*!
 *  @file TypeDefinition.hpp
 *  @brief Definition of fundamental types.
 */

#ifndef FVCODE3D_TYPEDEFNITION_HPP_
#define FVCODE3D_TYPEDEFNITION_HPP_

#include <FVCode3D/core/BasicType.hpp>
#include <FVCode3D/geometry/Point3D.hpp>
#include <FVCode3D/core/Chrono.hpp>
#include "GetPot"

namespace FVCode3D
{
  
  //! Type for a std::function<Real(Point3D)>
  typedef std::function<Real(Point3D)> Func;
  //! For time dependent problems
  typedef std::function<Real(Real,Point3D)> Funct;
  
} // namespace FVCode3D

#endif /* FVCODE3D_TYPEDEFNITION_HPP_ */
