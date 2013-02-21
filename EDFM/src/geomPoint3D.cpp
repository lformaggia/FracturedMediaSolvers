 /*!
 *	@file geomPoint3D.cpp
 *	@brief Struct for Point and Vector in 3D space (definition). 
 *
 *	@author Luca Turconi <lturconi@gmail.com>
 *  @date 22-09-2012
 *
 */

#include "geomPoint3D.hpp"

namespace Geometry
{
	
	// --------------------   Struct Point3D   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================
	Point3D::Point3D() : x(), y(), z() {}

	Point3D::Point3D(const Real & a, const Real & b, const Real & c) : x(a), y(b), z(c) {}

	Point3D::Point3D(const Point3D & p) : x(p.x), y(p.y), z(p.z) {}

	Point3D::~Point3D() {}
	
// ==================================================
// Methods
// ==================================================
	Point3D Point3D::cross(const Point3D & p) const
	{
		Point3D tmp( y*p.z-z*p.y ,
					 z*p.x-x*p.z ,
					 x*p.y-y*p.x );
		return tmp;
	}

	void Point3D::showMe(std::ostream  & out) const
	{
		out << "Type = Point3D/Vector3D : ";
		out << " ( " << x << " , "
				<< y << " , "
				<< z << " )" << std::endl;
	}

// ==================================================
// Operators
// ==================================================
	bool Point3D::operator==( const Point3D & b)
	{
		if (fabs(x-b.x)<1.0e-5 &&fabs(y-b.y)<1.0e-5 &&fabs(z-b.z)<1.0e-5 )		
		return true;
		else
		return false;
	}

	Point3D & Point3D::operator+=( const Point3D & b)
	{
		x += b.x;
		y += b.y;
		z += b.z;
		
		return *this;
	}
	
	Point3D & Point3D::operator-=( const Point3D & b)
	{
		x -= b.x;
		y -= b.y;
		z -= b.z;
		
		return *this;
	}
	
	Point3D & Point3D::operator*=( const Real & a)
	{
		x *= a;
		y *= a;
		z *= a;
		
		return *this;
	}
		
	Point3D & Point3D::operator/=( const Real & a)
	{
		x /= a;
		y /= a;
		z /= a;
		
		return *this;
	}
	
// ==================================================
// External Operators
// ==================================================
	Point3D const operator+(const Point3D & p1, const Point3D & p2)
	{
		Point3D sum(p1);
		sum += p2;
		return sum;
	}
	
	Point3D const operator-(const Point3D & p1, const Point3D & p2)
	{
		Point3D diff(p1);
		diff -= p2;
		return diff;
	}
	
	Point3D const operator*(const Real & a, const Point3D & p)
	{
		Point3D tmp(p);
		tmp *= a;
		return tmp;
	}
	
	Point3D const operator*(const Point3D & p, const Real & a)
	{
		Point3D tmp(p);
		tmp *= a;
		return tmp;
	}
	
	Point3D const operator/(const Point3D & p, const Real & a)
	{
		Point3D tmp(p);
		tmp /= a;
		return tmp;
	}
	
	std::ostream& operator<<(std::ostream & ostr, const Point3D & p)
	{
		ostr << " ( " << p.x << " , "
				<< p.y << " , "
				<< p.z << " )" << std::flush;
		return ostr;
	}
	
	bool operator<(const Point3D & A, const Point3D & B)
	{
		if( A.x-B.x < -eps )
			return 1;
		if( A.x-B.x > eps )
			return 0;
		if( std::fabs(A.x-B.x) < eps)
		{
			if( A.y-B.y < -eps )
				return 1;
			if( A.y-B.y > eps )
				return 0;
			if( std::fabs(A.y-B.y) < eps)
			{
				if( A.z-B.z < -eps )
					return 1;
				if( A.z-B.z > eps )
					return 0;
				if( std::fabs(A.z-B.z) < eps)
					return 0;
			}
		}
		return 0;
	}

	
	
} // namespace Geometry
