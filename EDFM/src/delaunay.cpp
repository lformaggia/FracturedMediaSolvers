#include "delaunay.hpp"
#include <algorithm>
#include <iterator>
#include <stack>
#include <limits>
#include <utility>
#include <set>
#include <cmath>
// unnamed namespace to hide globals used only in this compilation unit
namespace
{
  //!Indicates that an identifier is null.
  const unsigned int notAnId = std::numeric_limits<unsigned int>::max();

  //! A function used only locally
  /*!
    It scales coordinates
  */
  void scaleCoord (std::vector<Geometry::Point2D>& coor)
  {
    using namespace Geometry;
    double xmin[2];
    double xmax[2];

    xmin[0] = coor[0].x;
    xmin[1] = coor[0].y;
    xmax[0] = coor[0].x;
    xmax[1] = coor[0].y;
    for (std::vector<Geometry::Point2D>::iterator i = coor.begin() + 1; i < coor.end(); ++i)
    {
      xmin[0] = std::min (xmin[0], i->x);
      xmin[1] = std::min (xmin[1], i->y);
      xmax[0] = std::max (xmax[0], i->x);
      xmax[1] = std::max (xmax[1], i->y);
    }
    double delta[2];
    delta[0] = xmax[0] - xmin[0];
    delta[1] = xmax[1] - xmin[1];
    Point2D orig (xmin[0], xmin[1]);
    //scale points to avoid problems
    for (std::vector<Point2D>::iterator i = coor.begin(); i < coor.end(); ++i)
      i->scale (orig, delta);
  }// end scale

  //! A Helper functor used only locally
  class TakeOutElements
  {
  public:
    TakeOutElements (unsigned int test) : M_test (test) {};
    bool operator() (Geometry::IdTriplet const& i) const
    {
      return i[0] >= M_test || i[1] >= M_test || i[2] >= M_test;
    }
  private:
    unsigned int M_test;
  };
  // end of TakeOutElements

  //! Mesh edges are just a couple of ints.
  typedef std::pair<unsigned int, unsigned int> BareEdge;

  //! Ordering for edges.
  /*!
    @note Ignoring orientation.
   */
  struct OrderEdges
  {
    bool operator () (BareEdge const& a, BareEdge const& b) const
    {
      unsigned int max_a = std::max (a.first, a.second);
      unsigned int max_b = std::max (b.first, b.second);
      if (max_a != max_b) return max_a < max_b;
      return std::min (a.first, a.second) < std::min (b.first, b.second);
    }
  };
}// end of unnamed namespace


namespace Geometry
{

  void delaun (std::vector<Point2D> const& points, std::vector<IdTriplet>& elements)
  {
    std::stack<unsigned int> stack_element;
    unsigned int numPoints (points.size() );
    std::vector<Point2D> coor;
    coor.reserve (numPoints + 3);
    coor.insert (coor.end(), points.begin(), points.end() );
    /// Scaling coordinates, better work in [0, 1]
    scaleCoord (coor);
    // Additional points to form the first triangle covering domain
    Point2D pa (-100.0, -100.0);
    Point2D pb ( 100.0, -100.0);
    Point2D pc (-10.0 , 100.0);
    coor.push_back (pa);
    coor.push_back (pb);
    coor.push_back (pc);
    // The first element
    IdTriplet v;
    std::vector<IdTriplet> elem;
    elem.reserve (2 * numPoints + 1);
    elem.push_back (make_triplet (numPoints, numPoints + 1, numPoints + 2) );
    // element2element map. For each element the adjacent element across a side
    // side i is the side with local numbering i-(i+1)%3
    std::vector<IdTriplet> sid;
    sid.reserve (3 * numPoints);
    sid.push_back (make_triplet (notAnId, notAnId, notAnId) );
    IdTriplet e;
    IdTriplet s;
    IdTriplet vnew;
    unsigned int enew[2];
    for (unsigned int p = 0; p < numPoints; ++p)
    {
      Point2D const& xp (coor[p]);
      // Locate element containing xp
      unsigned int t = triloc (xp, coor, elem, sid);
      if (t == notAnId)
      {
        std::string message;
        message = std::string ("cannot locate element. Error in") +
                  std::string ("__FILE__ line: __LINE__");
        throw std::runtime_error (message);
      }
      // New element by splitting current element
      e = sid[t]; //(a b c)
      v = elem[t]; // v1 v2 v3
      // numtri +1 +2 +3
      enew[0] = elem.size();
      enew[1] = enew[0] + 1;
      // p v1 v2 -> replace t
      elem[t] = make_triplet (p, v[0], v[1]);
      // numtri+ 2 a numtri +1
      sid[t] = make_triplet (enew[1], e[0], enew[0]);
      // p v2 v3 -> new element
      elem.push_back (make_triplet (p, v[1], v[2]) );
      // t b  numtri + 2
      sid.push_back (make_triplet (t, e[1], enew[1]) );
      // p v3 v1 -> new element
      elem.push_back (make_triplet (p, v[2], v[0]) );
      // numtri + 1 c t
      sid.push_back (make_triplet (enew[0], e[2], t) );
      if (e[0] != notAnId)
      {
        stack_element.push (t);
      }
      if (e[1] != notAnId)
      {
        sid[e[1]][edg (e[1], t, sid)] = enew[0];
        stack_element.push (enew[0]);
      }
      if (e[2] != notAnId)
      {
        sid[e[2]][edg (e[2], t, sid)] = enew[1];
        stack_element.push (enew[1]);
      }
      unsigned int l, erl, era, erb;
      unsigned int r, a, b, c;
      unsigned int v1, v2, v3;
      Point2D x1, x2, x3;
      while (!stack_element.empty() )
      {
        l = stack_element.top();
        stack_element.pop();
        // get element adjacent to second side of (l)
        r = sid[l][1];
        //find corresponding side in other element (r)
        erl = edg (r, l, sid);
        era = (erl + 1) % 3;
        erb = (era + 1) % 3;
        // Get nodes of element r in the right order
        v1  = elem[r][erl];
        v2  = elem[r][era];
        v3  = elem[r][erb];
        x1  = coor[v1];
        x2  = coor[v2];
        x3  = coor[v3];
        // Test if swapping is good
        if (swap (x1, x2, x3, xp) )
        {
          a = sid[r][era];
          b = sid[r][erb];
          c = sid[l][2];
          elem[l][2] = v3;
          sid[l][1] = a;
          sid[l][2] = r;
          elem[r] = make_triplet (p, v3, v1);
          sid[r] = make_triplet (l, b, c);
          if (a != notAnId)
          {
            sid[a][edg (a, r, sid)] = l;
            stack_element.push (l);
          }
          if (b != notAnId)
          {
            stack_element.push (r);
          }
          if (c != notAnId)
          {
            sid[c][edg (c, l, sid)] = r;
          }
        } // end if side has swapped

      } // end while on stack_element
    } // end for Points

    // We know the number of elements by Euler's formula
    if (elem.size() != 2 * points.size() + 1)
    {
      throw std::runtime_error ("Error while constructing Delaunay mesh. __FILE__  __LINE__");
    }
    // Remove all elements whose point id is greater or equal the number of points
    TakeOutElements oper (points.size() );
    // remove_if does not actually remove.
    std::vector<IdTriplet>::iterator new_end = std::remove_if (elem.begin(), elem.end(), oper);
    // COpy the good elements into elements.
    elements.clear();
    // MAke sure we have space to avoid useless allocations
    elements.reserve (std::distance (elem.begin(), new_end) );
    // Copy
    std::copy (elem.begin(), new_end, std::back_inserter (elements) );
  } // end delaunay

  unsigned int edg (
    unsigned int const& element1,
    unsigned int const& element2,
    std::vector<IdTriplet> const& sides)
  {
    IdTriplet const& edge = sides[element1];
    for (unsigned int i = 0; i < 3; ++i)
    {
      if (edge[i] == element2) return i;
    }
    throw std::runtime_error ("Not finding side in edg");
  }// end edg

  unsigned int triloc (Point2D const& xp, std::vector<Point2D> const& coor,
                       std::vector<IdTriplet> const& elem,
                       std::vector<IdTriplet> const& sides)
  {
    unsigned int counter (0), v1 (0), v2 (0);
    double d1, d2, d3, d4;
    unsigned int t (elem.size() - 1);
    bool test (false);
    do
    {
      if (++counter > 2 * elem.size() )
        throw std::runtime_error ("Cannot find point in mesh. File __FILE__, Line __LINE__");
      for (unsigned int i = 0; i < 3; ++i)
      {
        // get element side
        v1 = elem[t][i      ];
        v2 = elem[t][ (i + 1) % 3];
        // coefficients of the line defined by the side
        d1 = coor[v1].y - xp.y;
        d2 = coor[v2].x - xp.x;
        d3 = coor[v2].y - xp.y;
        d4 = coor[v1].x - xp.x;
        // Check if point is outside!
        test = d1 * d2 > d3 * d4;
        if (test)
        {
          t = sides[t][i];
          break; // get out of the loop
        }
      }// end for  on i
    }
    while (test);
    return t;
  }// end triloc

  bool swap (Point2D const& x1, Point2D const& x2, Point2D const& x3, Point2D const& xp)
  {
    Point2D x13, x23, x1p, x2p;
    double cosa, cosb, sina, sinb;
    bool result;
    x13 = x1 - x3;
    x23 = x2 - x3;
    x1p = x1 - xp;
    x2p = x2 - xp;
    cosa = x13.x * x23.x + x13.y * x23.y;
    cosb = x1p.x * x2p.x + x1p.y * x2p.y;
    if ( (cosa >= 0.) && (cosb >= 0.0) )
      result = false;
    else
    {
      if ( (cosa < 0.0) && (cosb < 0) )
        result = true;
      else
      {
        sina = x13.x * x23.y - x23.x * x13.y;
        sinb = x2p.x * x1p.y - x1p.x * x2p.y;
        result = (sina * cosb + sinb * cosa) < 0.0;
      }
    }
    return result;
  }// end swap

  std::vector<unsigned int>
  orderedBoundaryPoints (std::vector<IdTriplet> const& elements)
  {
    OrderEdges ordering;
    std::set<BareEdge, OrderEdges> edges (ordering);
    std::set<BareEdge, OrderEdges>::iterator l;
    BareEdge e;
    for (unsigned int i = 0; i < elements.size(); ++i)
    {
      for (unsigned int j = 0; j < 3 ; ++j)
      {
        e.first  = elements[i][j      ];
        e.second = elements[i][ (j + 1) % 3];
        if ( (l = edges.find (e) ) == edges.end() )
        {
          edges.insert (e);
        }
        else
        {
          edges.erase (l);
        }
      }// end for j
    }// end for i
    // Now I have all boundary edges. I need to order them.
    std::vector<unsigned int> result;
    result.reserve (edges.size() );
    l = edges.begin();
    // Start inserting the second point of first edge
    unsigned int next = l->second;
    result.push_back (next);
    edges.erase (l);
    bool ok;
    // Not the nicest algorithm: it is O(n^2), n being the number of
    // boundary edges. But I do not know how to do better.
    while (!edges.empty() )
    {
      ok = false;
      for (l = edges.begin(); l != edges.end(); ++l)
      {
        // Add second point of the edge whose first point
        // is equal to next. Adjourn next, delete side.
        if (l->first == next)
        {
          next = l->second;
          result.push_back (next);
          edges.erase (l);
          ok = true;
          break;
        }// end if
      }// end for l
      if (!ok) throw std::runtime_error ("Wrong mesh, __FILE__");
    }// end while
    return result;
  }// end  orderedBoundaryPoints


  double area (Point2D const& a, Point2D const& b, Point2D const& c)
  {
    const double& p1 = a.x;
    const double& p2 = b.x;
    const double& p3 = c.x;
    const double& q1 = a.y;
    const double& q2 = b.y;
    const double& q3 = c.y;
    return p2 * q3 - p3 * q2 - p1 * q3 + p3 * q1 + p1 * q2 - p2 * q1;
  }// end area

  double diameter2 (Point2D const& a, Point2D const& b, Point2D const& c)
  {
    const double  p1 = a.x - c.x;
    const double  p2 = b.x - a.x;
    const double  p3 = c.x - b.x;
    const double  q1 = a.y - c.y;
    const double  q2 = b.y - a.y;
    const double  q3 = c.y - b.y;
    // in C++11 you may use max with arbitrary number of args.
    return std::max (std::max (p1 * p1 + q1 * q1, p2 * p2 + q2 * q2), p3 * p3 + q3 * q3);
  }// end area

  bool Aligned::operator() (Point2D const& a, Point2D const& b, Point2D const& c) const
  {
    return std::fabs (area (a, b, c) ) < M_tol * diameter2 (a, b, c);
  }

  std::vector<bool>
  Aligned::decimated_list (const std::vector<Point2D> points) const
  {
    std::vector<bool> eliminated (points.size(), false);
    Aligned const& self (*this);
    unsigned int i (0), j (1), k (2);
    for (unsigned int ii = 0; ii < points.size(); ++ii)
    {
      if (self (points[i], points[j], points[k]) )
      {
        eliminated[j] = true;
        j = k;
      }
      else
      {
        i = j;
        j = k;
      }
      ++k;
      k = k % points.size();
    }// end for
    return eliminated;
  }// end decimated_list

  std::vector<Point2D>
  Aligned::decimate (const std::vector<Point2D> points) const
  {
    std::vector<bool> eliminated = decimated_list (points);
    unsigned int count (0);
    for (unsigned int i = 0; i < points.size(); ++i)
      if (!eliminated[i]) ++count;
    std::vector<Point2D> result;
    if (count == 0) return result;
    result.reserve (count);
    for (unsigned int i = 0; i < points.size(); ++i)
      if (!eliminated[i]) result.push_back (points[i]);
    return result;
  }// end decimate

} // end namespace Geometry
