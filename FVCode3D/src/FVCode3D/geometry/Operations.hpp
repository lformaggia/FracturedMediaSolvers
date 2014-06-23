/*!
 *  @file operations.hpp
 *  @brief Some useful operations.
 */

#ifndef OPERATIONS_HPP_
#define OPERATIONS_HPP_

#include <FVCode3D/core/TypeDefinition.hpp>

namespace FVCode3D
{

class Point3D;

//! Compute the normal moving from three points that lie on a plane (use the right-hand rule)
/*!
 * @param A first point
 * @param B second point
 * @param C third point
 * @return signed normal
 */
Point3D computeNormal(const Point3D & A, const Point3D & B, const Point3D & C);

//! Rotate a point of an angle
/*!
 * @param p point to rotate
 * @param angleDeg angle of rotation (in Degrees)
 * @return the rotated point
 */
Point3D rotateOf(const Point3D & p, const Real angleDeg);

//! Compute the area of a triangle
/*!
 * @param A first point of the triangle
 * @param B second point of the triangle
 * @param C third point of the triangle
 * @return unsigned area of the triangle
 */
Real triangleArea(const Point3D & A, const Point3D & B, const Point3D & C);

//! Compute the volume of a tetrahedron
/*!
 * @param nodes vector of points that define the tetrahedron
 * @return unsigned volume of the tetrahedron
 */
Real tetrahedronVolume(const std::vector<Point3D> & nodes);

//! Functor for computing the minimum between two types
template <typename T>
struct less
{
	bool operator() (const T & x, const T & y) const { return x < y; }
};

//! Functor for computing the minimum between two generic std::pair
template <typename T>
struct less< std::pair<T,T> >
{
	bool operator() (const std::pair<T,T> & x, const std::pair<T,T> & y) const
	{
		if( x.first < y.first )
			return true;
		if( x.first == y.first )
		{
			if( x.second < y.second )
				return true;
			else
				return false;
		}
		else
			return false;
	}
};

}// namespace FVCode3D

#endif /* OPERATIONS_HPP_ */


