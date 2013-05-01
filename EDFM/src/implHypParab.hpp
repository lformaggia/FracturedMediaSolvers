#ifndef HH_IMPLHYPPARAB_HH
#define HH_IMPLHYPPARAB_HH
#include "geomPoint3D.hpp"
/*! \file implHypParab.hpp 

  We implement an implicit representation of a bilinear surface patch.
*/
namespace Geometry{
  /*! A simple struct that holds the 4 points defining the bilinear patch.

    The points are ordered in consecutive order along the patch boundary and
    such that the normal will be oriented by the right hand rule.
   */
    
  struct bilinerarPatchPoints{
    Point3D P00;
    Point3D P10;
    Point3D P11;
    Point3D P01;
  };
 
 /*! Implicit representation of a bilinear patch.
    
    Returns zero if (x,y,z) is on the patch. The returned value is greater
    than zero if the point is in the positive sempispace according to the 
    patch normal direction. 

    P (u,v) = (1-u) (1-v) P(0,0) + u (1-v) P(1,0) + v(1-u) P(0,1) + u v P(1,1)
    P=P00+(P10-P00)u+(P01-P00)v+(P11+P00-P10-P01)uv 
    P= [a0,b0,c0]+[a1,b1,c1]u + [a2,b2,c2]v + [a3,b3,c3]uv 
    
   
  */
    
  double implHypParab(const bilinerarPatchPoints & patch, double const x, double const y, double const z);
}
#endif
