/*!
 * @file Point3D.cpp
 * @brief Class for a 3D Point (definitions).
 */

#include "Point3D.hpp"
#include "CoordinateSystem.hpp"

namespace Geometry{

void Point3D::normalize(){
	Real n=norm();
	M_x = M_x/n;
	M_y = M_y/n;
	M_z = M_z/n;
}

void Point3D::setValues(const Real x, const Real y, const Real z){
	M_x = x;
	M_y = y;
	M_z = z;
}

void Point3D::linearTransform(const Real scaling, const Real xShift, const Real yShift, const Real zShift){

	M_x = M_x * scaling + xShift;
	M_y = M_y * scaling + yShift;
	M_z = M_z * scaling + zShift;

}

Point3D Point3D::convertInLocalCoordinate(const CoordinateSystem3D & coordSys, const Point3D & origin) const{
	Point3D temp;
	temp = *this - origin;
	Real xloc( dotProduct( temp , coordSys.getU() ) );
	Real yloc( dotProduct( temp , coordSys.getV() ) );
	Real wloc( dotProduct( temp , coordSys.getW() ) );

	temp.setValues(xloc,yloc,wloc);
	return temp;
}

Point3D & Point3D::operator=(const Point3D & p){
	if ( this != &p ){
		M_x = p.M_x;
		M_y = p.M_y;
		M_z = p.M_z;
	}
	return *this;
}

Point3D & Point3D::operator+=(const Point3D & p){
	M_x += p.M_x;
	M_y += p.M_y;
	M_z += p.M_z;
	return *this;
}

Point3D & Point3D::operator-=(const Point3D & p){
	M_x -= p.M_x;
	M_y -= p.M_y;
	M_z -= p.M_z;
	return *this;
}

Point3D & Point3D::operator*=(const Real & r){
	M_x *= r;
	M_y *= r;
	M_z *= r;
	return *this;
}

Point3D & Point3D::operator/=(const Real & r){
	M_x /= r;
	M_y /= r;
	M_z /= r;
	return *this;
}

Real Point3D::operator[](const UInt coord) const{

	switch(coord){

	case 0:
		return M_x;
		break;
	case 1:
		return M_y;
		break;
	case 2:
		return M_z;
		break;
	default:
		// TODO throw exception
		return M_x;

	}

}

Real & Point3D::operator[](const UInt coord){

	switch(coord){

	case 0:
		return M_x;
		break;
	case 1:
		return M_y;
		break;
	case 2:
		return M_z;
		break;
	default:
		// TODO throw exception
		return M_x;

	}

}

Point3D operator+(const Point3D & p1, const Point3D & p2){
	Point3D sum(p1);
	sum += p2;
	return sum;
}

Point3D operator-(const Point3D & p1, const Point3D & p2){
	Point3D diff(p1);
	diff -= p2;
	return diff;
}

Point3D operator-(const Point3D & p)
{
	return p*(-1);
}

Real operator*(const Point3D & p1, const Point3D & p2)
{
	return dotProduct(p1,p2);
}

Point3D operator*(const Point3D & p, const Real r)
{
	Point3D mul(p);
	mul *= r;
	return mul;
}

Point3D operator/(const Point3D & p, const Real r){
	Point3D div(p);
	div /= r;
	return div;
}

Real dotProduct(const Point3D & p1, const Point3D & p2){
	Real p(p1.M_x * p2.M_x + p1.M_y * p2.M_y + p1.M_z * p2.M_z);
	return p;
}

Real innerAngleRad(const Point3D & p1, const Point3D & p2){
	Real theta = dotProduct(p1,p2) / (p1.norm() * p2.norm());
	theta = theta > 1 ? 0 : acos(theta);
	return theta > _PI_ ? 2 * _PI_ - theta : theta;
}

Real innerAngleDeg(const Point3D & p1, const Point3D & p2){
	Real theta = dotProduct(p1,p2) / (p1.norm() * p2.norm());
	theta = theta > 1 ? 0 : acos(theta) / _PI_ * 180;
	return theta > 180 ? 360 - theta : theta;
}

Real distance(const Point3D & p1, const Point3D & p2){
	return (p1-p2).norm();
}

Point3D crossProduct(const Point3D & p1, const Point3D & p2){
	Point3D p(p1.M_y*p2.M_z - p1.M_z*p2.M_y,
			p1.M_z*p2.M_x - p1.M_x*p2.M_z,
			p1.M_x*p2.M_y - p1.M_y*p2.M_x);
	return p;
}

std::ostream & operator<<(std::ostream & os, const Point3D & p){
	os<<"( "<<p.M_x<<" , "<<p.M_y<<" , "<<p.M_z<<" )";
	return os;
}

}//namespace Geometry