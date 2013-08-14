#ifndef HH__DELAUNAY_HH
#define HH__DELAUNAY_HH
#include <cstring>
#include <vector>
#include <limits>
#include <stdexcept>
#include "geomPoint2D.hpp"
#include "tolerances.hpp"
/*!
  @file delaunay.hpp
  @author Luca Formaggia
  
  Utilities for triangulating sets of 2D points and
  and post-processing triangulations.
 */
namespace Geometry
{
  //! Holds a triplet of values.
  template<typename T>
  struct Triplet
  {
    T & operator[](unsigned int i) {return M_Id[i]; }
    T operator[](unsigned int i) const {return M_Id[i]; }
    Triplet() {};
    Triplet(const Triplet<T> & rhs)
    {
      std::memcpy(this->M_Id,rhs.M_Id,3*sizeof(T));
    };
    Triplet & operator =(const Triplet<T> & rhs)
    {
      if (this != &rhs)
	std::memcpy(this->M_Id,rhs.M_Id,3*sizeof(T));
      return *this;
    }
  private:
    T M_Id[3];
  };
  //! A helper function to make a triplet.
  template <typename T>
  inline Triplet<T> make_triplet(T a, T b, T c)
  {
    Triplet<T> tmp;
    tmp[0]=a;
    tmp[1]=b;
    tmp[2]=c;
    return tmp;
  }

  //! A triplet of unsigned int
  typedef Triplet<unsigned int> IdTriplet;

  //! A very basic Delaunay triangulator.
  /*!  
    It implements the algorithm by Lawson. Not the most effective
    one on average, so this routine should be used for small to
    moderate (<=10^4) number of points.
    
    \parameter points A vector of two dimensional points.
    \parameter elements The vector of element connectivity.
    \throw runtime_error
  */
  void delaun(std::vector<Point2D> const & points, std::vector<IdTriplet> & elements);
  
  //! Finds which edge is shared between two  given elements.
  /*!
    The algoritm looks in sides[element1][] to find element2. The algorithm throw a runtime_error if the element is not found.

    @parameter element1 The first element. Is the one whose edgeto element table is looked up.
    @parameter element2 The second element to be matched.
    @parameter sides    The element2element map.
    @return The found element side. 
    @throw  runtime_error.
   */
  unsigned int edg(unsigned int const & element1, unsigned int const & element2, std::vector<IdTriplet> const & sides);

  //! Locates a point in a mesh.
  /*!  Given a delaunay mesh triangulating a convex region and a
    point, the function locates the element containing the point.
    
    @parameter xp. Point to be located;
    @parameter coor. Vector with point coordinates.
    @parameter elem. Mesh connectivity. As a vector of triplets.
    @parameter sides. The element2element map. For each element the adjacent elements.
   */
  unsigned int triloc(Point2D const & xp, std::vector<Point2D> const & coor, std::vector<IdTriplet> const & elem, std::vector<IdTriplet> const & sides);
  
  //! Test if swapping sides is good.
  /*!
    The configuration 1 2-3 4 is compared with 1 4-2 3 by flipping a side.
    The function return true if  the second one is Delaunay satisfying.
   */
  bool swap(Point2D const & x1, Point2D const & x2, Point2D const & x3, Point2D const & x4);
  //! Extracts points on the boundary of triangulation.
  /*!
    Given a triangulation it returns the list of boundary points,
    well ordered.
    
    \param elements. The connectivity matrix of a valid 2D
    triangulation of a simply connected domain.
    
    \return the connectivity The ids of the points at the boundary, ordered
    
    \pre \c elements should contain the connectivity matrix of a valid 2D
    triangulation of a simply connected domain. That is, its boundary
    must be formed by a \b single \b polygonal \b line. Furthermore all mesh
    elements must be ordeed consistently. The ordering of the points
    in output depends on the ordering of the mesh element.
  */
  std::vector<unsigned int> 
  orderedBoundaryPoints(std::vector<IdTriplet> const & elements);

  //! Area (with sign) of a triangle defined by 3 points.
  /*!
    Sign is positive if points a,b and c are counterclockwise oriented.
   */
  double area(Point2D const & a, Point2D const & b, Point2D const & c);

  //! Square of the diameter of a triangle defined by 3 points.
  double diameter2(Point2D const & a, Point2D const & b, Point2D const & c);
  
  //! Detects if points are aligned
  struct Aligned{
    //! Default contructor. Tolerance may be varied.
    Aligned(Real const tol=EDFM_Tolerances::ALIGNMENT_TOLERANCE):M_tol(tol){}
    //! Returns true if the points are on a line.
    /*
      \pre Points cannot be all coincident. This condition is a precondition, so it is not checked by the
      method.
     */
    bool operator()(Point2D const & a, Point2D const & b, Point2D const & c) const;
    //! Decimates a vector of 2D points.
    /*!
      It identifies the points that are aligned canceling
      the middle one.
      
      \param point The points to be decimated.
      \return vector of bool of the same lenght as the input. true if point is decimated.
    */
    std::vector<bool> decimated_list(const std::vector<Point2D> points) const;
    //! Decimates a vector of 2D points.
    /*!
      It elimiminates points that are aligned canceling
      the middle one.
      
      \param point The points to be decimated.
      \return The decimated points.
    */
    std::vector<Point2D> decimate(const std::vector<Point2D> points) const;

  private:
    Real const M_tol;
  };

}// end namespace Geometry.
#endif
