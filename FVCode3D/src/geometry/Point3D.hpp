/*!
 * @file Point3D.hpp
 * @brief Class for a 3D Point.
 */

#ifndef POINT3D_HPP_
#define POINT3D_HPP_

#include "core/basic_type.hpp"

namespace Geometry{

class CoordinateSystem3D;

//! Class that defines a 3D Point and some operations.
/*!
 * @class Point3D
 * This class defines a point as three Real coordinates.
 * It also implements some operation between two points or between a point and a scalar value.
 */
class Point3D{
public:

	Point3D(): M_x(0.) ,M_y(0.) , M_z(0.) {};

	Point3D(Real x, Real y, Real z): M_x(x) ,M_y(y) , M_z(z) {};

	Point3D(const Point3D & p): M_x(p.x()), M_y(p.y()), M_z(p.z()) {};

	Real x() const { return M_x; };

	Real y() const { return M_y; };

	Real z() const { return M_z; };

	Real norm() const { return std::sqrt( M_x*M_x + M_y*M_y + M_z*M_z );}

	void normalize();

	void setValues(const Real x, const Real y, const Real z);

	void linearTransform(const Real scaling, const Real xShift, const Real yShift, const Real zShift);

	Point3D convertInLocalCoordinate(const CoordinateSystem3D & coordSys, const Point3D & origin) const;

	Point3D & operator=(const Point3D & p);

	Point3D & operator+=(const Point3D & p);

	Point3D & operator-=(const Point3D & p);

	Point3D & operator*=(const Real & r);

	Point3D & operator/=(const Real & r);

	Real operator[](const UInt coord) const;

	Real & operator[](const UInt coord);

	static const Real & getTolerance()
	    { return Point3D::S_tolerance; }

	static void setTolerance(const Real tolerance)
	    { Point3D::S_tolerance = tolerance; }

	friend Point3D operator+(const Point3D & p1, const Point3D & p2);

	friend Point3D operator-(const Point3D & p1, const Point3D & p2);

	friend Point3D operator-(const Point3D & p);

	friend Real operator*(const Point3D & p1, const Point3D & p2);

	friend Point3D operator*(const Point3D & p, const Real r);

	friend Point3D operator*(const Real r, const Point3D & p);

	friend Point3D operator/(const Point3D & p, const Real r);

	friend Real dotProduct(const Point3D & p1, const Point3D & p2);

	friend Real innerAngleRad(const Point3D & p1, const Point3D & p2);

	friend Real innerAngleDeg(const Point3D & p1, const Point3D & p2);

	friend Real distance(const Point3D & p1, const Point3D & p2);

	friend Point3D crossProduct(const Point3D & p1, const Point3D & p2);

	friend std::ostream & operator<<(std::ostream & os, const Point3D & p);

private:

	//! x-coordinate
	Real M_x;
	//! y-coordinate
	Real M_y;
	//! z-coordinate
	Real M_z;

	//! Relative tolerance
	static Real S_tolerance;

};

bool operator<(const Point3D & p1, const Point3D & p2);

}//namespace Geometry

#endif /* POINT3D_HPP_ */
