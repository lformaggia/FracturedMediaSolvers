 /*!
 *	@file Typedefinition3D.hpp
 *	@brief Class for unstructured mesh.
 *
 *	@author Francesco Della Porta 
 *
 */ 


#ifndef TYPEDEFNITION3D_HPP_
#define TYPEDEFNITION3D_HPP_

#include<CGAL/Point_3.h>
#include<CGAL/Vector_3.h>
#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>

namespace Geometry
{

//! CGAL::kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

//! Type for real numbers of CGAL::kernel
typedef Kernel::FT Real;

//! Type for integer numbers
typedef long unsigned int UInt;
	
//! Type for 3D vector
typedef CGAL::Vector_3< Kernel > Vector3D;

//! Type for 3D points
typedef CGAL::Point_3< Kernel > Point3D;

//! Type for 3D segments
typedef CGAL::Segment_3< Kernel > Segment3D;

//! Type for Mesh3D
class Mesh3D
{
public:
	class Facet3D;
	class Cell3D;
};

//! Type for FractureNetwork3D
class FractureNetwork3D
{};

} // namespace Geometry

#endif /* TYPEDEFNITION_HPP_ */
