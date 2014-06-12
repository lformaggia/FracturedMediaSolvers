/*!
 * @file CoordinateSystem.cpp
 * @brief Class for a 3D Cartesian coordinate system (definitions).
 */

#include "CoordinateSystem.hpp"

namespace Geometry
{

void CoordinateSystem3D::computeCartesianCoordinateSystem(const Point3D & z)
{
	M_w = z;
	M_w.normalize();

	//find the maximum of M_w components
	UInt max = 0;
	Point3D temp1;
	Point3D temp2;
	Point3D x2;
	Point3D x3;

	if ( std::abs(M_w.y()) > std::abs(M_w.x()) )
		if ( std::abs(M_w.z()) > std::abs(M_w.y()) )
			max = 2;
		else
			max = 1;
	else if ( std::abs(M_w.z()) > std::abs(M_w.x()) )
		max = 2;

	if ( max==0 )
	{
		x2.setValues(0., 1., 0.);
		x3.setValues(0., 0., 1.);
	}
	else if ( max==1 )
	{
		x2.setValues(1., 0., 0.);
		x3.setValues(0., 0., 1.);
	}
	else if ( max==2 )
	{
		x2.setValues(1., 0., 0.);
		x3.setValues(0., 1., 0.);
	}

	// u assignement
	temp1 = M_w;
	temp1.linearTransform( dotProduct(x2,M_w), 0. , 0. , 0. );
	M_u = x2 - temp1;
	M_u.normalize();

	// v assignement
	temp1 = M_w;
	temp2 = M_u;

	temp1.linearTransform( dotProduct(x3,M_w) , 0. , 0. , 0. );
	temp2.linearTransform( dotProduct(x3,M_u) , 0. , 0. , 0. );

	M_v = x3 - temp1;
	M_v -= temp2;

	M_v.normalize();
}

}// namespace Geometry
