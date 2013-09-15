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
#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Different Kernel
//#include<CGAL/Cartesian.h>
//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>

//! Type for real numbers
typedef double Real;

//! Type for real numbers
typedef unsigned int UInt;

namespace Geometry
{
	
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

//typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
//typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt Kernel;
	
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
