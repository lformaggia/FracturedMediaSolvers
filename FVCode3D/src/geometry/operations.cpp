/*!
 *  @file operations.cpp
 *  @brief Some useful operations (definitions).
 */

#include "operations.hpp"
#include "geometry/Point3D.hpp"

namespace FVCode3D{

Point3D computeNormal( const Point3D & A, const Point3D & B, const Point3D & C)
{
	Point3D n( crossProduct(A-C,B-C) );
	n.normalize();
	return n;
}

Point3D rotateOf(const Point3D & p, const Real angleDeg)
{
	Real angleRad = angleDeg / 180 * _PI_;
	Point3D r1(std::cos(angleRad), std::sin(angleRad), 0.);
	Point3D r2(-std::sin(angleRad), std::cos(angleRad), 0.);

	Real x, y;
	x = dotProduct(r1,p);
	y = dotProduct(r2,p);

	return Point3D(x,y,p.z());
}

Real triangleArea(const Point3D & A, const Point3D & B, const Point3D & C)
{
	return crossProduct(B-A,C-A).norm() / 2.;
}

Real tetrahedronVolume(const std::vector<Point3D> & nodes)
{
	return std::fabs( dotProduct( nodes[2]-nodes[3] , crossProduct(nodes[0]-nodes[3],nodes[1]-nodes[3]) ) ) / 6.;
}

}// namespace FVCode3D
