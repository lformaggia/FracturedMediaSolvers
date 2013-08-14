 /*!
 *	@file geomPoint2D.cpp
 *	@brief Struct for Point and Vector in 2D space (definition). 
 *
 *	@author Luca Turconi <lturconi@gmail.com>
 *  @date 22-09-2012
 *
 */

#include "geomPoint2D.hpp"
#include "tolerances.hpp"

namespace Geometry
{
	
	// --------------------   Struct Point2D   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================
	Point2D::Point2D() : x(0.), y(0.) {}

	Point2D::Point2D(const Real & a, const Real & b) : x(a), y(b){}
  
  Point2D:: Point2D(std::vector<Real> const &v): x(v[0]),y(v[1]) {};
        

	Point2D::Point2D(const Point2D & p) : x(p.x), y(p.y){}

	Point2D::~Point2D() {}
	
// ==================================================
// Methods
// ==================================================
	double Point2D::cross(const Point2D & p) const
	{
	  return x*p.y-y*p.x;
	}

	void Point2D::showMe(std::ostream  & out) const
	{
		out << "Type = Point2D/Vector2D : ";
		out << " ( " << x << " , "
				<< y << " )" << std::endl;
	}

  void Point2D::scale(Point2D const & orig, double factors[])
  {
    x=(x-orig.x)/factors[0];
    y=(y-orig.y)/factors[1];
  }

  void Point2D::scaleBack(Point2D const & orig, double factors[])
  {
    x=factors[0]*x+orig.x;
    y=factors[1]*y+orig.y;
  }
// ==================================================
// Operators
// ==================================================
  bool Point2D::operator==( const Point2D & b)const
  {
    if (fabs(x-b.x)==0. &&fabs(y-b.y)==0.)		
      return true;
    else
      return false;
  }
  
  Point2D & Point2D::operator+=( const Point2D & b)
  {
    x += b.x;
    y += b.y;
    
    return *this;
  }
  
  Point2D & Point2D::operator-=( const Point2D & b)
  {
    x -= b.x;
    y -= b.y;
    
    return *this;
  }
  
  Point2D & Point2D::operator*=( const Real & a)
  {
    x *= a;
    y *= a;
    
    return *this;
  }
  
  Point2D & Point2D::operator/=( const Real & a)
  {
    x /= a;
    y /= a;
    
    return *this;
  }
	
// ==================================================
// External Operators
// ==================================================
  Point2D const operator+(const Point2D & p1, const Point2D & p2)
  {
    Point2D sum(p1);
    sum += p2;
    return sum;
  }
  
  Point2D const operator-(const Point2D & p1, const Point2D & p2)
  {
    Point2D diff(p1);
    diff -= p2;
    return diff;
  }
  
  Point2D const operator*(const Real & a, const Point2D & p)
  {
    Point2D tmp(p);
    tmp *= a;
    return tmp;
  }
  
  Point2D const operator*(const Point2D & p, const Real & a)
  {
    Point2D tmp(p);
    tmp *= a;
    return tmp;
  }
  
  Point2D const operator/(const Point2D & p, const Real & a)
  {
    Point2D tmp(p);
    tmp /= a;
    return tmp;
  }
  
  std::ostream& operator<<(std::ostream & ostr, const Point2D & p)
  {
    ostr << " ( " << p.x << " , "
	 << p.y << " )" << std::flush;
    return ostr;
  }
  
  bool operator<(const Point2D & A, const Point2D & B)
  {
    const Real eps=0.0;
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
      }
    return 0;
  }
  
} // namespace Geometry
