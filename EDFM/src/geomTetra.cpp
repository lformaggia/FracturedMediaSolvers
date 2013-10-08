/*!
*  @file geomTetra.cpp
*  @brief Class for Tetra in 3D space (definition).
*
*
*/

#include<cmath>
#include<limits>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<vector>
#include<algorithm>
#include<cassert>
#include<cstdlib>
#include "geomTetra.hpp"

namespace Geometry
{

  // --------------------   Class Tetra  --------------------

  // ==================================================
  // Constructors & Destructor
  // ==================================================
  Tetra::Tetra() : M_pA(), M_pB(), M_pC(), M_pD() {}

  Tetra::Tetra (const Point3D& a, const Point3D& b, const Point3D& c, const Point3D& d) :
    M_pA (a), M_pB (b), M_pC (c), M_pD (d)
  {
    ;
  }

  Tetra::Tetra (const Tetra& t) : M_pA (t.A() ), M_pB (t.B() ), M_pC (t.C() ), M_pD (t.D() )
  {
    ;
  }

  Tetra::~Tetra() {}

  Point3D Tetra::getPoint (const gmm::size_type& quale)
  {
    //LF: better safe than sorry
    assert (quale >= 0 && quale <= 3);
    if (quale == 0) return M_pA;
    if (quale == 1) return M_pB;
    if (quale == 2) return M_pC;
    if (quale == 3) return M_pD;
    std::cerr << " Errore in getPoint" << std::endl;
    std::exit (1);

  }

  Triangle Tetra::getFace (const Point3D& p1, const Point3D& p2 , const Point3D& p3)
  {
    Triangle t (p1, p2, p3);
    return t;
  }

  Triangle Tetra::getFace (const gmm::size_type& i1, const gmm::size_type& i2, const gmm::size_type& i3)
  {
    Triangle t (this->getPoint (i1), this->getPoint (i2), this->getPoint (i3) );
    return t;
  }

  std::vector<Point3D> Tetra::getGaussNodes()
  {
    Point3D p;
    std::vector<Point3D> nodigauss;
    p = 0.585410196624969 * M_pA + 0.138196601125011 * M_pB + 0.138196601125011 * M_pC + 0.138196601125011 * M_pD;
    nodigauss.push_back (p);
    p = 0.585410196624969 * M_pB + 0.138196601125011 * M_pA + 0.138196601125011 * M_pC + 0.138196601125011 * M_pD;
    nodigauss.push_back (p);
    p = 0.585410196624969 * M_pC + 0.138196601125011 * M_pB + 0.138196601125011 * M_pA + 0.138196601125011 * M_pD;
    nodigauss.push_back (p);
    p = 0.585410196624969 * M_pD + 0.138196601125011 * M_pB + 0.138196601125011 * M_pC + 0.138196601125011 * M_pA;
    nodigauss.push_back (p);
    return nodigauss;
  }

  std::vector<Real> Tetra::getGaussWeights()
  {
    std::vector<Real>  pesigauss (4, 0.25);
    return pesigauss;
  }


} // namespace Geometry
