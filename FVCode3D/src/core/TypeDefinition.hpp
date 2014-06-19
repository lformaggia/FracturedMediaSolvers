/*!
 *  @file TypeDefinition.hpp
 *	@brief Definition of fundamental types.
 */

#ifndef TYPEDEFNITION_HPP_
#define TYPEDEFNITION_HPP_

#include "core/BasicType.hpp"
#include "geometry/Point3D.hpp"
#include "core/Chrono.hpp"
#include "GetPot"

namespace FVCode3D
{

typedef	std::function<Real(Point3D)> Func;

} // namespace FVCode3D

#endif /* TYPEDEFNITION_HPP_ */
