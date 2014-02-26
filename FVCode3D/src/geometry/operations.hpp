/*!
 *  @file operations.hpp
 *  @brief Some useful operations.
 */

#ifndef OPERATIONS_HPP_
#define OPERATIONS_HPP_

#include "core/TypeDefinition.hpp"

namespace Geometry{

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

	//! Compute the volume of a tetrahedron
	/*!
	 * @param nodes vector of points that define the tetrahedron
	 * @return unsigned volume of the tetrahedron
	 */
	Real tetrahedronVolume(const std::vector<Point3D> & nodes);

}// namespace Geometry

#endif /* OPERATIONS_HPP_ */


