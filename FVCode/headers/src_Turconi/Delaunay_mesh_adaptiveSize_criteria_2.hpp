/*!
 *	@file Delaunay_mesh_adaptiveSize_criteria_2.hpp
 *	@brief Class for unstructured triangular mesh in 2D space.
 *
 *	@author Luca Turconi <luca.turconi@moxoff.com>
 *
 */

// Copyright (c) 2003-2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/releases/CGAL-4.1-branch/Mesh_2/include/CGAL/Delaunay_mesh_size_criteria_2.h $
// $Id: Delaunay_mesh_size_criteria_2.h 70936 2012-08-01 13:29:16Z lrineau $
// 
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_DELAUNAY_MESH_ADAPTIVESIZE_CRITERIA_2_HPP
#define CGAL_DELAUNAY_MESH_ADAPTIVESIZE_CRITERIA_2_HPP

#include "TypeDefinition.hpp"
#include "Fracture2D.hpp"

#include <CGAL/Mesh_2/Face_badness.h>
#include <CGAL/Delaunay_mesh_criteria_2.h>

#include <utility>
#include <ostream>
#include <cmath>

namespace CGAL {

template <class CDT>
class Delaunay_mesh_adaptiveSize_criteria_2 : 
    public virtual Delaunay_mesh_criteria_2<CDT>
{
protected:
  typedef typename CDT::Geom_traits Geom_traits;
  Geometry::FractureNetwork2D * M_fnPtr;
  Real M_hmin;
  Real M_hmax;
  Real M_alpha;
  Real M_transitionRegion;
public:
  typedef Delaunay_mesh_criteria_2<CDT> Base;

  Delaunay_mesh_adaptiveSize_criteria_2() {}

  Delaunay_mesh_adaptiveSize_criteria_2(Geometry::FractureNetwork2D * const fn,
	  				const Real aspect_bound = 0.125, 
        	                        const Real hmin = 0,
        	                        const Real hmax = 0,
        	                        const Real transitionRegion = 0,
        	                        const Real alpha = 1,
        	                        const Geom_traits& traits = Geom_traits())
    : Base(aspect_bound, traits), M_fnPtr(fn), M_hmin(hmin), M_hmax(hmax), M_alpha(alpha) 
    {
    	if(transitionRegion>0)
    		M_transitionRegion = transitionRegion;
    	else
    		M_transitionRegion = std::fabs(hmax-hmin);
    }

  inline
  Real alpha() const { return M_alpha; }
  
  inline
  Real hmin() const { return M_hmin; }

  inline
  Real hmax() const { return M_hmax; }
  
  inline
  Real transitionRegion() const { return M_transitionRegion; }
  
  inline
  Geometry::FractureNetwork2D * fnPtr() const { return M_fnPtr; }

  inline
  void set_alpha(const Real sb) { M_alpha = sb; }

  inline
  void set_hmin(const Real sb) { M_hmin = sb; }
  
  inline
  void set_hmax(const Real sb) { M_hmax = sb; }
  
  inline
  void set_transitionRegion(const Real sb) { M_transitionRegion = sb; }
  
  inline
  void set_fnPtr(const Geometry::FractureNetwork2D * fn) { M_fnPtr = fn; }

  // first: squared_minimum_sine
  // second: size
  struct Quality : public std::pair<Real, Real>
  {
    typedef std::pair<Real, Real> Base;

    Quality() : Base() {};
    Quality(Real _sine, Real _size) : Base(_sine, _size) {}

    const Real & size() const { return second; }
    const Real & sine() const { return first; }

    // q1<q2 means q1 is prioritised over q2
    // ( q1 == *this, q2 == q )
    bool operator<(const Quality & q) const
    {
      if( size() > 1 )
	if( q.size() > 1 )
	  return ( size() > q.size() );
	else
	  return true; // *this is big but not q
      else
	if( q.size() >  1 )
	  return false; // q is big but not *this
      return( sine() < q.sine() );
    }

    std::ostream& operator<<(std::ostream& out) const
    {
      return out << "(size=" << size()
		 << ", sine=" << sine() << ")";
    }
  };

  class Is_bad: public Base::Is_bad
  {
  protected:
//    Real squared_size_bound; // squared size bound on edge length
    Geometry::FractureNetwork2D * const M_fnPtr;
    const Real M_hmin;
    const Real M_hmax;
    const Real M_alpha;
    const Real M_transitionRegion;
  public:
    typedef typename Base::Is_bad::Point_2 Point_2;

    Is_bad(Geometry::FractureNetwork2D * const fn,
	   const Real aspect_bound, 
           const Real hmin,
           const Real hmax,
           const Real transitionRegion,
           const Real alpha,
           const Geom_traits& traits)
      : Base::Is_bad(aspect_bound, traits), M_fnPtr(fn),
        M_hmin(hmin), M_hmax(hmax), M_alpha(alpha), M_transitionRegion(transitionRegion) {}

    Mesh_2::Face_badness operator()(const Quality q) const
    {
      if( q.size() > 1 )
	return Mesh_2::IMPERATIVELY_BAD;
      if( q.sine() < this->B )
	return Mesh_2::BAD;
      else
	return Mesh_2::NOT_BAD;
    }

    Mesh_2::Face_badness operator()(const typename CDT::Face_handle& fh,
				    Quality& q) const
    {
      typedef typename CDT::Geom_traits Geom_traits;
      typedef typename Geom_traits::Compute_area_2 Compute_area_2;
      typedef typename Geom_traits::Compute_squared_distance_2
	Compute_squared_distance_2;

      Geom_traits traits; /** @warning traits with data!! */

      Compute_squared_distance_2 squared_distance = 
	traits.compute_squared_distance_2_object();

      const Point_2& pa = fh->vertex(0)->point();
      const Point_2& pb = fh->vertex(1)->point();
      const Point_2& pc = fh->vertex(2)->point();

      Real
	a = CGAL::to_double(squared_distance(pb, pc)),
	b = CGAL::to_double(squared_distance(pc, pa)),
	c = CGAL::to_double(squared_distance(pa, pb));
      
      Geometry::Point2D Ma( (pb.x()+pc.x())/2, (pb.y()+pc.y())/2 );
      Geometry::Point2D Mb( (pa.x()+pc.x())/2, (pa.y()+pc.y())/2 );
      Geometry::Point2D Mc( (pa.x()+pb.x())/2, (pa.y()+pb.y())/2 );
      
      Real max_sq_length; // squared max edge length
      Real second_max_sq_length;
      
      if(a<b)
	{
	  if(b<c) {
	    max_sq_length = c;
	    second_max_sq_length = b;
	  }
	  else { // c<=b
	    max_sq_length = b;
	    second_max_sq_length = ( a < c ? c : a );
	  }
	}
      else // b<=a
	{
	  if(a<c) {
	    max_sq_length = c;
	    second_max_sq_length = a;
	  }
	  else { // c<=a
	    max_sq_length = a;
	    second_max_sq_length = ( b < c ? c : b );
	  }
	}

      q.second = 0;
      
//      std::cout << " M_hmin = " << M_hmin << std::endl;
      
      if( M_hmin != 0 )
        {
	  Real aBound( M_hmin + (M_hmax - M_hmin) *
	  		std::pow( M_fnPtr->distance(Ma)/M_transitionRegion, M_alpha) );
	  Real bBound( M_hmin + (M_hmax - M_hmin) *
      			std::pow( M_fnPtr->distance(Mb)/M_transitionRegion, M_alpha) );
	  Real cBound( M_hmin + (M_hmax - M_hmin) *
	  		std::pow( M_fnPtr->distance(Mc)/M_transitionRegion, M_alpha) );
      
	  aBound *= aBound;
	  bBound *= bBound;
	  cBound *= cBound;
      
	  Real aQuality( a / aBound );
	  Real bQuality( b / bBound );
	  Real cQuality( c / cBound );
      	  
      	  if( aQuality>bQuality )
      	  {
      	  	if( aQuality>cQuality )
      	  		q.second = aQuality;
      	  	else
      	  		q.second = cQuality;
      	  }
      	  else // aQuality < bQuality
      	  {
      	  	if( bQuality>cQuality )
      	  		q.second = bQuality;
      	  	else
      	  		q.second = cQuality;
      	  }
      	  
	    // normalized by size bound to deal
	    // with size field
	  if( q.size() > 1 )
	    {
	      q.first = 1; // (do not compute sine)
	      return Mesh_2::IMPERATIVELY_BAD;
	    }
	}

      Compute_area_2 area_2 = traits.compute_area_2_object();

      double area = 2*CGAL::to_double(area_2(pa, pb, pc));

      q.first = (area * area) / (max_sq_length * second_max_sq_length); // (sine)
      
      if( q.sine() < this->B )
	return Mesh_2::BAD;
      else
	return Mesh_2::NOT_BAD;	
    }
  };

  Is_bad is_bad_object() const
  { return Is_bad(this->fnPtr(), this->bound(), this->hmin(), this->hmax(),
  		  this->transitionRegion(), this->alpha(),
                  this->traits /* from the bad class */); }
                  
};

} // end namespace CGAL

#endif
