 /*!
 *	@file geomBilinearSurface.hpp
 *	@brief Class for Bilinear Surface in 3D space.
 *
 *	@author Luca Turconi <lturconi@gmail.com>
 *  @date 22-09-2012
 *
 */
 
#ifndef GEOMBILINEARSURFACE_HPP_
#define GEOMBILINEARSURFACE_HPP_

#include<iostream>
#include<string>
#include<vector>

#include "TypeDefinition.hpp"
#include "geomPoint3D.hpp"
#include "geomSegment.hpp"
#include "geomTriangle.hpp"
#include "gmm/gmm.h"
#include "tolerances.hpp"
namespace Geometry
{
  /*!
    @class BilinearSurface
    	
    @author Luca Turconi <lturconi@gmail.com>
    	
    This class implements the concept of bilinear surface.
    It stores the four points generating the surface: M_pA, M_pB, M_pC, M_pD.
    The surface is created considering M_pA and M_pC as opposite points.
    	
    The class provides methods to evaluate the parametrization of the surface
    and its first and second derivatives, for different parameters values.
    	
    It is possible to compute the maximal curvature in each point of the surface,
    and to extract the maximal value of curvature on the whole surface.
		
    A tool of exact and approximated methods to compute the interesection between the surface
    and a segment are available.
		
    The method showMe allows to print the class informations.
    It is also possible to print the attributes of the surface by using << operator.
    	
    For 3D visualization use the exportVtk method.
    	
  */
  class BilinearSurface
  {
  public:
    //! @name Constructor & Destructor
    //@{
		
    //! Empty constructor
    BilinearSurface();
	
    //! Constructor, getting the extremal points
    /*!
     * Attention: when creating the surface, a and c are considered as opposite points.
     * @param a The first point
     * @param b The second point
     * @param c The third point
     * @param d The fourth point
     */
    BilinearSurface(const Point3D & a, const Point3D & b, const Point3D & c, const Point3D & d);
	
    //! Copy constructor
    /*!
     * @param b The bilinear surface copied in the new object
     */
    BilinearSurface(const BilinearSurface & b);
	
    //! Destructor
    virtual ~BilinearSurface();
	
    //@}
	
    //! @name Get Methods
    //@{
		
    //! Get point A
    /*!
     * @return The extreme point A
     */
    inline Point3D A() const { return M_pA; }
	
    //! Get point B
    /*!
     * @return The extreme point B
     */
    inline Point3D B() const { return M_pB; }
	
    //! Get point C
    /*!
     * @return The extreme point C
     */
    inline Point3D C() const { return M_pC; }
	
    //! Get point D
    /*!
     * @return The extreme point D
     */
    inline Point3D D() const { return M_pD; }
	
    //! Get surface maximum dimension
    /*!
     * @return The M_Lmax value
     */
    inline Real Lmax() const { return M_Lmax; }
	
    //! Get triangulation size
    /*!
     * @return The M_TriangulatedSurface size
     */
    inline UInt getTriangulationSize() const
    { return M_TriangulatedSurface.size(); }
	
    //@}
	
    //! @name Set Methods
    //@{
		
    //! Set point A
    /*!
     * @param p The new value for point A
     */
    inline void setA(const Point3D & p) { M_pA=p; }
	
    //! Set point B
    /*!
     * @param p The new value for point B
     */
    inline void setB(const Point3D & p) { M_pB=p; }	
	
    //! Set point C
    /*!
     * @param p The new value for point C
     */
    inline void setC(const Point3D & p) { M_pC=p; }	
	
    //! Set point D
    /*!
     * @param p The new value for point D
     */
    inline void setD(const Point3D & p) { M_pD=p; }	

    void rescaleZ(Real scale, Real shift);
	
    //! Set Lmax
    /*!
     * Set the maximum length of the surface
     */
    void setLmax();
	
    //@}
	
    //! @name Methods
    //@{
	
    //! Parametric equation
    /*!
     * It is the parametrization of the bilinear surface.
     * @param u The first parameter value (range [0,1])
     * @param v The second parameter value (range [0,1])
     * @return The point corresponding to the parameters value
     */
    inline Point3D param(const Real & u, const Real & v) const
    { return (1-u)*( (1-v)*M_pA+v*M_pD ) + u*( (1-v)*M_pB+v*M_pC ); }

    //! Inverse parametric equation
    /*!
     * Version that uses the normal.
     * It is just an approximation. It computes the inverse map by intersecting with 
     * a segment build from given point using the surface normal.
     *
     * @param p The 3D point on the surface of which we want to find the parametric coordinates.
     * @param normal the normal at the center of the surface
     * @return The parameter values.
     */

    std::vector<Real> inv_param(const Point3D & p ,const Point3D & normal) const;

    //! Inverse parametric equation
    /*!
     * Version that does not  use the normal.
     * It uses projection. Indeed it is just a wrapper to that method, using default parameters for the
     * Newton algorithm. You may just want to use the method projection directly to change tolerances.
     *
     * @param p The 3D point on the surface of which we want to find the parametric coordinates.
     * @return The parameter values.
     */
    std::vector<Real> inv_param(const Point3D & p) const;

    void makePlanar();

    Real areaFault();
    //! Parametric equation: first order partial derivative with respect to u
    /*!
     * It is the parametrization of the first derivative with respect to u
     * of the bilinear surface.
     * @param u The first parameter value (range [0,1])
     * @param v The second parameter value (range [0,1])
     * @return The point corresponding to the parameters value
     */
    inline Point3D S_u(const Real & u, const Real & v) const
    { return (1-v)*( M_pB-M_pA ) + v*( M_pC-M_pD ); }
	
    //! Parametric equation: first order partial derivative with respect to v
    /*!
     * It is the parametrization of the first derivative with respect to v
     * of the bilinear surface.
     * @param u The first parameter value (range [0,1])
     * @param v The second parameter value (range [0,1])
     * @return The point corresponding to the parameters value
     */
    inline Point3D S_v(const Real & u, const Real & v=0) const
    { return (1-u)*( M_pD-M_pA ) + u*( M_pC-M_pB ); }

    //! Normal vector
    /*!
     * It computes the normalized vector normal to the surface
     * at point associated with (u,v).
     * @return The normalized normal vector in (u,v)
     */
    inline Point3D normal(const Real & u, const Real & v) const
    { return S_u(u,v).cross(S_v(u,v)) / (S_u(u,v).cross(S_v(u,v)).norm()); }
	
    //! Parametric equation: mixed second order partial derivative
    /*!
     * It is the parametrization of the mixed second order partial derivative
     * of the bilinear surface.
     * @param u The first parameter value (range [0,1])
     * @param v The second parameter value (range [0,1])
     * @return The point corresponding to the parameters value
     */
    inline Point3D S_uv(const Real & u=0, const Real & v=0) const
    { return M_pA+M_pC-M_pB-M_pD; }
		
    //! Maximal curvature in (u,v)
    /*!
     * It computes the maximal curvature at point (u,v).
     * @param u The first parameter value (range [0,1])
     * @param v The second parameter value (range [0,1])
     * @return The maximal curvature at (u,v)
     */
    Real maxCurvature(const Real & u, const Real & v) const;
	
    //! Maximal curvature of the surface
    /*!
     * It computes the maximal curvature on the whole surface.
     * @return The maximal curvature on the whole surface
     */
    Real maxCurvature() const;
	
    //! Build surface Triangulation
    /*!
     * It creates a triangulation in order to approximate the surface
     * with an error less than the tolerance.
     * The triangulation are stored in M_TriangulatedSurface.
     * @param toll The approximation tolerance
     */
    void buildTriangulation(const Real & toll=eps);
	
    //! Erase surface Triangulation
    /*!
     * It erases the surface triangulation.
     */
    inline void cleanTriangulation()
    {
      M_Triangulated=0;
      M_TriangulatedSurface.clear();
    }
	
    //! Approximated test for intersection with a segment
    /*!
     * The surface is approximated with two triangles.
     * The test for intersection is performed between the segment and this triangles.
     * In the case of coplanar segment, the method will return FALSE.
     * (A single intersection does not exist!)
     * @param s The segmet to be tested
     * @param stdDivision A bool parameter to chose the way to split the surface.
     * @return TRUE if the segment has an intersection with a triangle
     * 		FALSE if the intersection does not exist
     */	
    bool approxIsIntersectedBy(const Segment & s, const bool & stdDivision=1) const;
	
    //! Compute exact intersection with a segment
    /*!
     * It computes the intersection between the surface and a given segment.
     * The nonlinear problem is solved using the standard Newton's method.
     * If the intersection does not exist, this method will return the point (NaN,NaN,NaN).
     * Equally, in the case of coplanar segment, the method will return (NaN,NaN,NaN).
     * (A single intersection does not exist!)
     * @param s The segment
     * @param toll The Newton's method tolerance
     * @param maxIter The Newton's method iteration limit
     * @return The intersection point or (NaN,NaN,NaN) if the intersection does not exist.
     */	
    Point3D newtonIntersectionWith(const Segment & s, const Real & toll=eps, const UInt & maxIter=60) const;
    //! Version that returns also the point in the parameter plane.
    void newtonIntersectionWith_uv(const Segment & s, std::vector<Real> &, const Real & toll=eps, const UInt & maxIter=60) const;
    
    //! Compute approximated intersection with a segment
    /*!
     * The surface is approximated with two triangles.
     * It computes the intersection between the surface approximation and a given segment.
     * If the intersection does not exist, this method will return the point (NaN,NaN,NaN).
     * Equally, in the case of coplanar segment, the method will return (NaN,NaN,NaN).
     * (A single intersection does not exist!)
     * @param s The segment,
     * @param stdDivision A bool parameter to chose the way to split the surface.
     * @return The intersection point or (NaN,NaN,NaN) if the intersection does not exist.
     */	
    Point3D approxIntersectionWith(const Segment & s, const bool & stdDivision=1) const;
	
    //! Approximated test with exact computation of intersection with a segment
    /*!
     * The surface is approximated with two triangles.
     * The method performs a test for the existence of intersections with the segment
     * using the surface approximation.
     * If the test return TRUE, it will compute the exact intersection between the surface
     * and a given segment.
     * The nonlinear problem is solved using the standard Newton's method.
     * If the intersection does not exist, this method will return the point (NaN,NaN,NaN).
     * Equally, in the case of coplanar segment, the method will return (NaN,NaN,NaN).
     * (A single intersection does not exist!)
     * @param s The segment
     * @param toll The Newton's method tolerance
     * @param maxIter The Newton's method iteration limit
     * @return The intersection point or (NaN,NaN,NaN) if the intersection does not exist.
     */	
    Point3D approxNewtonIntersectionWith(const Segment & s,
					 const Real & toll=eps, const UInt & maxIter=60) const;
	
    //! Compute intersection between a given segment and the surface triangulation
    /*!
     * It computes the intersection between the surface triangulation and a given segment.
     * If the intersection does not exist, this method will return the point (NaN,NaN,NaN).
     * Equally, in the case of coplanar segment, the method will return (NaN,NaN,NaN).
     * (A single intersection does not exist!)
     * @param s The segment
     * @param toll The tolerance of the method
     * @return The intersection point or (NaN,NaN,NaN) if the intersection does not exist.
     */	
    Point3D approxRefinedIntersectionWith(const Segment & s, const Real & toll=eps);
	
    //! Export in vtk format
    /*!
     * It generates the vtk file, for the 3D visualization of the bilinear surface.
     * (Use paraview to open the vtk file)
     * @param filename The name of the vtk file created by this method
     * @return TRUE -> operation ended correctly
     FALSE -> an error occurred
    */	
    //! Project a 3D point onto the surface.
    /*!
      It uses Newton algorithm to project a point onto the surface.
      No backtracking anabled so you should make sure that the point
      is not too far from the surface.
      \param s  The point to be projected
      \param uv The coordinates in the parametric plane of the projected point.
      \param tol The tolerance.
      \param maxIter Max number of iterations.
      \return The projected 3D point.
     */
    Point3D projection(const Point3D & s,std::vector<Real> & uv, 
		       const Real & toll=eps,
		       const UInt & maxIter=60) const;

    bool isIntersectedBy(BilinearSurface  altra, std::vector<Point3D> &) const; 

    bool isEdgeIntersectedBy(gmm::size_type quale, Segment S, Point3D &) const;

    bool exportVtk(const std::string & filename) const;
	
    //! Export triangulation in vtk format
    /*!
     * It generates the vtk file, for the 3D visualization of the surface triangulation.
     * (Use paraview to open the vtk file)
     * @param filename The name of the vtk file created by this method
     * @return TRUE -> operation ended correctly
     FALSE -> an error occurred
    */	
    bool exportTriangulationVtk(const std::string & filename) const;
	
    //! Display general information about the content of the class
    /*!
     * List of things displayed in the class
     * @param out Specify the output format (std::cout by default)
     */
    void showMe(std::ostream & out=std::cout) const;
	
    //@}
	
  private:
    Point3D M_pA;
    Point3D M_pB;
    Point3D M_pC;
    Point3D M_pD;
	
    Real M_Lmax;
    bool M_Triangulated;
    Real M_TriangulationToll;
    std::vector<Triangle> M_TriangulatedSurface;
  };

  //! @name External Operators
  //@{
	
  //! The insertion operator
  /*!
   * @param ostr The stream object on which the action is performed.
   * 			This is the first parameter of the global functions,
   * 			and represents the object to the left of the operator,
   * 			i.e. the object on which the extraction operation is performed.
   * @param b The bilinear surface inserted on the stream.
   * @return The stream object on which the action is performed (ostr)
   */
  std::ostream& operator<<(std::ostream & ostr, const BilinearSurface & b);
	
  //@}
	
} // namespace Geometry

#endif /* GEOMBILINEARSURFACE_HPP_ */
