 /*!
 *  @file TypeDefnition.hpp
 *	@brief Definition of fundamental scalar and geometrical types.
 *
 *	@author Luca Turconi <luca.turconi@moxoff.com>
 *
 */

#ifndef TYPEDEFNITION_HPP_
#define TYPEDEFNITION_HPP_

#include<CGAL/Point_2.h>
#include<CGAL/Segment_2.h>
#include<CGAL/Triangle_2.h>
#include<CGAL/Vector_2.h>
#include<CGAL/Line_2.h>
#include<CGAL/Polygon_2.h>

// Different Kernels
#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>


//! Type for real numbers and CGAL::kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT Real;

//! Type for integer numbers
typedef long unsigned int UInt;

namespace Geometry
{	
//! Type for 2D points
typedef CGAL::Point_2< Kernel > Point2D;

//! Type for 2D segments
typedef CGAL::Segment_2< Kernel > Segment2D;

//! Type for 2D segments
typedef CGAL::Triangle_2< Kernel > Triangle2D;

//! Type for 2D vectors
typedef CGAL::Vector_2< Kernel > Vector2D;

//! Type for 2D lines
typedef CGAL::Line_2< Kernel > Line2D;

//! Type for 2D polygon
typedef CGAL::Polygon_2< Kernel > Polygon2D;

} // namespace Geometry

#endif /* TYPEDEFNITION_HPP_ */
