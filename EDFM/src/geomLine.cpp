/*!
*  @file geomLine.cpp
*  @brief Class for Line in 3D space (definition).
*
*  @author Luca Turconi <lturconi@gmail.com>
*  @date 22-09-2012
*
*/

#include "geomLine.hpp"

namespace Geometry
{

  // --------------------   Class Line   --------------------

  // ==================================================
  // Constructors & Destructor
  // ==================================================
  Line::Line() : M_pA(), M_pB() {}

  Line::Line (const Real& xA, const Real& yA, const Real& zA,
              const Real& xB, const Real& yB, const Real& zB) :
    M_pA (xA, yA, zA), M_pB (xB, yB, zB) {}

  Line::Line (const Point3D& a, const Point3D& b) : M_pA (a), M_pB (b) {}

  Line::Line (const Segment& s) : M_pA ( s.A() ), M_pB ( s.B() ) {}

  Line::Line (const Line& l) : M_pA ( l.A() ), M_pB ( l.B() ) {}

  Line::~Line() {}

  // ==================================================
  // Methods
  // ==================================================
  Point3D Line::getPointAtX (const Real& x) const
  {
    Vector3D v (M_pB - M_pA);
    Real alpha = (x - M_pA.x) / v.x;
    if ( v.x == 0 && x != M_pA.x )
    {
      std::cerr << " *** Error: line has no points with x=" << x << " *** " << std::endl;
      return M_pA;
    }
    Point3D p (x, alpha * v.y + M_pA.y, alpha * v.z + M_pA.z);
    return p;
  }

  Point3D Line::getPointAtY (const Real& y) const
  {
    Vector3D v (M_pB - M_pA);
    Real alpha = (y - M_pA.y) / v.y;
    if ( v.y == 0 && y != M_pA.y )
    {
      std::cerr << " *** Error: line has no points with y=" << y << " *** " << std::endl;
      return M_pA;
    }
    Point3D p (alpha * v.x + M_pA.x, y, alpha * v.z + M_pA.z);
    return p;
  }

  Point3D Line::getPointAtZ (const Real& z) const
  {
    Vector3D v (M_pB - M_pA);
    Real alpha = (z - M_pA.z) / v.z;
    if ( v.z == 0 && z != M_pA.z )
    {
      std::cerr << " *** Error: line has no points with z=" << z << " *** " << std::endl;
      return M_pA;
    }
    Point3D p (alpha * v.x + M_pA.x, alpha * v.y + M_pA.y, z);
    return p;
  }

  void Line::showMe (std::ostream&   out) const
  {
    out << "Type = Line : " << std::endl;
    out << " M_pA : ";
    M_pA.showMe();
    out << " M_pB : ";
    M_pB.showMe();
  }

  // ==================================================
  // Operators
  // ==================================================
  std::ostream& operator<< (std::ostream& ostr, const Line& l)
  {
    ostr << " [ " << l.A() << " , "
         << l.B() << " ]" << std::flush;
    return ostr;
  }

} // namespace Geometry