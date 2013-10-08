/*!
*  @file geomSegment.cpp
*  @brief Class for Line in 3D space (definition).
*
*  @author Luca Turconi <lturconi@gmail.com>
*  @date 22-09-2012
*
*/

#include<limits>
#include<fstream>
#include <algorithm>
#include <iomanip>
#include <vector>
#include "geomSegment.hpp"
#include "geomTriangle.hpp"

namespace Geometry
{

  // --------------------   Class Segment   --------------------

  // ==================================================
  // Constructors & Destructor
  // ==================================================
  Segment::Segment() : M_pA(), M_pB() {}

  Segment::Segment (const Point3D& a, const Point3D& b) : M_pA (a), M_pB (b) {}

  Segment::Segment (const Segment& s) : M_pA (s.A() ), M_pB (s.B() ) {}

  Segment::~Segment() {}

  // ==================================================
  // Methods
  // ==================================================

  bool Segment::intersectTheSegment (const Segment& S, Point3D& risultato) const
  {
    Point3D d (M_pB - M_pA);
    Point3D m (S.A() - S.B() );
    Point3D P (S.A() - M_pA);
    Real t, s;
    /*std::cout << "----------------"<<std::endl;

    std::cout << M_pA<<std::endl;
    std::cout << M_pB<<std::endl;
    std::cout << S.A()<<std::endl;
    std::cout << S.B()<<std::endl;
    */
    if (d.x * m.y - m.x * d.y != 0)
    {
      t = (P.x * m.y - m.x * P.y) / (d.x * m.y - m.x * d.y);
      s = (d.x * P.y - P.x * d.y) / (d.x * m.y - m.x * d.y);
      if (t >= 0 && t <= 1 && s >= 0 && s <= 1 && fabs (d.z * t + m.z * s - P.z) < 1.0e-5 * d.norm() )
      {
        risultato = d * t + M_pA;
        //std::cout << risultato<<std::endl;
        return true;
      }
    }
    if (d.x * m.z - m.x * d.z != 0)
    {
      t = (P.x * m.z - m.x * P.z) / (d.x * m.z - m.x * d.z);
      s = (d.x * P.z - P.x * d.z) / (d.x * m.z - m.x * d.z);
      if (t >= 0 && t <= 1 && s >= 0 && s <= 1 && fabs (d.y * t + m.y * s - P.y) < 1.0e-5 * d.norm() )
      {
        risultato = d * t + M_pA;
        //std::cout << risultato<<std::endl;
        return true;
      }
    }
    if (d.y * m.z - m.y * d.z != 0)
    {
      t = (P.y * m.z - m.y * P.z) / (d.y * m.z - m.y * d.z);
      s = (d.y * P.z - P.y * d.z) / (d.y * m.z - m.y * d.z);
      if (t >= 0 && t <= 1 && s >= 0 && s <= 1 && fabs (d.x * t + m.x * s - P.x) < 1.0e-5 * d.norm() )
      {
        risultato = d * t + M_pA;
        //std::cout << risultato<<std::endl;
        return true;
      }
    }
    return false;




  }

  bool Segment::intersectThePlaneOf (const Triangle& t) const
  {
    Vector3D v1 ( M_pA - t.A() );
    Vector3D v2 ( M_pB - t.A() );

    Vector3D v3 ( this->param (0.5) - t.A() );

    if ( std::fabs ( v1.dot (t.normal() ) ) <= eps * t.Lmax()  &&
         std::fabs ( v2.dot (t.normal() ) ) <= eps * t.Lmax()  &&
         std::fabs ( v3.dot (t.normal() ) ) <= eps * t.Lmax()  )
    {
      // coplanar segment
      return 0;
    }

    if ( std::fabs ( v1.dot (t.normal() ) ) <= eps * t.Lmax() )
      return 1;

    if ( std::fabs ( v2.dot (t.normal() ) ) <= eps * t.Lmax() )
      return 1;

    if ( v1.dot (t.normal() ) *v2.dot (t.normal() )
         <= eps * t.Lmax() )
      return 1;

    return 0;
  }

  Point3D Segment::intersectionWithThePlaneOf (const Triangle& t) const
  {
    Point3D inter (std::numeric_limits<Real>::quiet_NaN(),
                   std::numeric_limits<Real>::quiet_NaN(),
                   std::numeric_limits<Real>::quiet_NaN() );

    Vector3D v1 ( M_pA - t.A() );
    Vector3D v2 ( M_pB - t.A() );

    Vector3D v3 ( this->param (0.5) - t.A() );

    if ( std::fabs ( v1.dot (t.normal() ) ) <= eps * t.Lmax() &&
         std::fabs ( v2.dot (t.normal() ) ) <= eps * t.Lmax() &&
         std::fabs ( v3.dot (t.normal() ) ) <= eps * t.Lmax() )
    {
      // coplanar segment
      return inter;
    }

    if ( std::fabs ( v1.dot (t.normal() ) ) <= eps * t.Lmax() )
      // M_pA is the intersection
      inter = M_pA;

    if ( std::fabs ( v2.dot (t.normal() ) ) <= eps * t.Lmax() )
      // M_pB is the intersection
      inter = M_pB;

    if ( v1.dot (t.normal() ) *v2.dot (t.normal() )
         <= eps * t.Lmax() )
    {
      // compute the intersection...
      Real coeff = std::fabs (v1.dot (t.normal() ) );
      if ( v1.dot (t.normal() ) > 0 )
        coeff /= (M_pA - M_pB).dot (t.normal() );
      if ( v1.dot (t.normal() ) < 0 )
        coeff /= (M_pB - M_pA).dot (t.normal() );

      inter = M_pA + (M_pB - M_pA) * coeff;
    }

    return inter;
  }

  Real  Segment::inv_param (const Point3D& p) const
  {
    Point3D vAP (p.x - M_pA.x, p.y - M_pA.y, p.z - M_pA.z);
    Point3D vAB (M_pB.x - M_pA.x, M_pB.y - M_pA.y, M_pB.z - M_pA.z);
    return vAP.dot (vAB) / (vAB.x * vAB.x + vAB.y * vAB.y + vAB.z * vAB.z);

  }


  bool Segment::isIn (Point3D const& P) const
  {
    Real t1, t2, t3;
    Real const& toll (EDFM_Tolerances::POINT_IN_SEGMENT);
    t1 = (fabs (M_pA.x - M_pB.x) > toll * fabs (M_pA.x) ) ? (P.x - M_pB.x) / (M_pA.x - M_pB.x) : 0.5;
    t2 = (fabs (M_pA.y - M_pB.y) > toll * fabs (M_pA.y) ) ? (P.y - M_pB.y) / (M_pA.y - M_pB.y) : 0.5;
    t3 = (fabs (M_pA.z - M_pB.z) > toll * fabs (M_pA.z) ) ? (P.z - M_pB.z) / (M_pA.z - M_pB.z) : 0.5;
    //
    if (t1 >= 0 && t1 <= 1 && t2 >= 0 && t2 <= 1 && t3 >= 0 && t3 <= 1 )
    {
      return true;
    }
    else
    {
      return false;
    }

  }

  bool Segment::exportVtk (const std::string& fileName) const
  {
    std::fstream filestr;

    filestr.open (fileName.c_str(), std::ios_base::out);

    if (filestr.is_open() )
    {
      std::cout << std::endl << " File: " << fileName << ", successfully opened";
    }
    else
    {
      std::cerr << std::endl << " *** Error: file not opened *** " << std::endl << std::endl;
      return  0;
    }

    std::cout << std::endl << " Exporting segment in Vtk format... " << std::endl;

    UInt nCells = 1;
    UInt nPoints = 2;
    UInt CellType = 3; // for VTK_HEXAHEDRON

    // Header
    filestr << "# vtk DataFile Version 3.1" << std::endl;
    filestr << "this is a file created for Paraview" << std::endl;
    filestr << "ASCII" << std::endl;
    filestr << "DATASET UNSTRUCTURED_GRID" << std::endl;
    filestr << std::endl; // The fifth line is empty.

    filestr << std::scientific << std::setprecision (10);

    // Pointdata
    filestr << "POINTS " << nPoints << " double" << std::endl;
    filestr << M_pA.x << " " << M_pA.y
            << " " << M_pA.z << std::endl;
    filestr << M_pB.x << " " << M_pB.y
            << " " << M_pB.z << std::endl;
    filestr << std::endl;

    // Celldata
    filestr << "CELLS " << nCells << " " << 3 << std::endl;
    filestr << "2 0 1" << std::endl;
    filestr << std::endl;

    filestr << "CELL_TYPES " << nCells << std::endl;
    filestr << CellType << std::endl;
    filestr << std::endl;

    filestr.close();

    return 1;
  }

  void Segment::showMe (std::ostream&   out) const
  {
    out << "Type = Segment : " << std::endl;
    out << " M_pA : ";
    M_pA.showMe();
    out << " M_pB : ";
    M_pB.showMe();
  }

  // ==================================================
  // Operators
  // ==================================================
  std::ostream& operator<< (std::ostream& ostr, const Segment& s)
  {
    ostr << " [ " << s.A() << " , "
         << s.B() << " ]" << std::flush;
    return ostr;
  }

  bool operator< (const Segment& s1, const Segment& s2)
  {
    // s1: Segment Orientation
    Point3D s1Inf, s1Sup;
    if (s1.B() < s1.A() )
    {
      s1Inf = s1.B();
      s1Sup = s1.A();
    }
    else
    {
      s1Inf = s1.A();
      s1Sup = s1.B();
    }
    // s2: Segment Orientation
    Point3D s2Inf, s2Sup;
    if (s2.B() < s2.A() )
    {
      s2Inf = s2.B();
      s2Sup = s2.A();
    }
    else
    {
      s2Inf = s2.A();
      s2Sup = s2.B();
    }

    // operator< implementation
    if ( ! ( s1Inf < s2Inf ) && ! ( s2Inf < s1Inf ) )
      return s1Sup < s2Sup;
    return s1Inf < s2Inf;
  }


  // --------------------   Class Segment2D   --------------------


  // ==================================================
  // Constructors & Destructor
  // ==================================================
  Segment2D::Segment2D() : M_pA(), M_pB() {}

  Segment2D::Segment2D (const Point2D& a, const Point2D& b) : M_pA (a), M_pB (b) {}

  Segment2D::Segment2D (const Segment2D& s) : M_pA (s.A() ), M_pB (s.B() ) {}

  Segment2D::~Segment2D() {}

  // ==================================================
  // Methods
  // ==================================================

  bool Segment2D::intersectTheSegment (const Segment2D& S,
                                       Point2D& risultato) const
  {
    Point2D S1_pA (this->M_pA);
    Point2D S1_pB (this->M_pB);
    Point2D S2_pA (S.A() );
    Point2D S2_pB (S.B() );
    Real bbox[4];
    // Find bounding box
    bbox[0] = std::min (S1_pA.x, (std::min (S1_pB.x, std::min (S2_pA.x, S2_pB.x) ) ) );
    bbox[1] = std::min (S1_pA.y, (std::min (S1_pB.y, std::min (S2_pA.y, S2_pB.y) ) ) );
    bbox[2] = std::max (S1_pA.x, (std::max (S1_pB.x, std::max (S2_pA.x, S2_pB.x) ) ) );
    bbox[3] = std::max (S1_pA.y, (std::max (S1_pB.y, std::max (S2_pA.y, S2_pB.y) ) ) );
    // Origin and scaling factors
    Point2D origin (bbox[0], bbox[1]);
    Real factors[2];
    factors[0] = std::max (EDFM_Tolerances::SEGMENT2D_INTERSECTION_TOLERANCE, bbox[2] - bbox[0]);
    factors[1] = std::max (EDFM_Tolerances::SEGMENT2D_INTERSECTION_TOLERANCE, bbox[3] - bbox[1]);
    // Change points so that we are operating in the unitary box
    S1_pA.scale (origin, factors);
    S1_pB.scale (origin, factors);
    S2_pA.scale (origin, factors);
    S2_pB.scale (origin, factors);

    // Build matrix
    Real A[2][2];
    A[0][0] = S1_pB.x - S1_pA.x;
    A[0][1] = S2_pA.x - S2_pB.x;
    A[1][0] = S1_pB.y - S1_pA.y;
    A[1][1] = S2_pA.y - S2_pB.y;
    Real det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
    // Parallel segments.
    if (std::fabs (det) < EDFM_Tolerances::SEGMENT2D_INTERSECTION_TOLERANCE)
    {
      Point2D normal (S1_pB.y - S1_pA.y, S1_pA.x - S1_pB.x);
      normal /= normal.norm();
      Real dist = std::fabs (Point2D (S1_pA - S2_pA).dot (normal) );
      if (dist > EDFM_Tolerances::POINT_IN_SEGMENT)
      {
        return false;
      }
      // I need to check all points!
      if (this->isIn (S.A() ) )
      {
        risultato = S.A();
        return true;
      }
      if (this->isIn (S.B() ) )
      {
        risultato = S.B();
        return true;
      }
      if (S.isIn (this->A() ) )
      {
        risultato = this->A();
        return true;
      }
      if (S.isIn (this->B() ) )
      {
        risultato = this->B();
        return true;
      }
      return false;
    }
    // Not parallel segments.
    double t[2];
    t[0] = ( A[1][1] * (S2_pA.x - S1_pA.x) - A[0][1] * (S2_pA.y - S1_pA.y) ) / det;
    t[1] = (-A[1][0] * (S2_pA.x - S1_pA.x) + A[0][0] * (S2_pA.y - S1_pA.y) ) / det;
    //if parametric coordinates outside [0,1] no intersection
    if (
      t[0] < -EDFM_Tolerances::POINT_IN_SEGMENT   ||
      t[0] > 1 + EDFM_Tolerances::POINT_IN_SEGMENT ||
      t[1] < -EDFM_Tolerances::POINT_IN_SEGMENT   ||
      t[1] > 1 + EDFM_Tolerances::POINT_IN_SEGMENT)
    {
      return false;
    }
    // It should be irrelevant how we compute the point
    risultato = t[0] * S1_pB + (1.0 - t[0]) * S1_pA;
    // Scale back to physical coordinates.
    risultato.scaleBack (origin, factors);
    return true;
  }


  Real  Segment2D::inv_param (const Point2D& p) const
  {
    Point2D vAP (p.x - M_pA.x, p.y - M_pA.y);
    Point2D vAB (M_pB.x - M_pA.x, M_pB.y - M_pA.y);
    return vAP.dot (vAB) / (vAB.x * vAB.x + vAB.y * vAB.y);
  }

  bool Segment2D::isIn (Point2D const& P) const
  {
    Real t = this->inv_param (P);
    if (
      t < -EDFM_Tolerances::POINT_IN_SEGMENT   ||
      t > 1 + EDFM_Tolerances::POINT_IN_SEGMENT) return false;
    Point2D Check = this->param (t);
    Real dist = Segment2D (P, Check).length();
    if (dist < EDFM_Tolerances::POINT_IN_SEGMENT * this->length() )
      return true;
    else
      return false;
  }



} // namespace Geometry
