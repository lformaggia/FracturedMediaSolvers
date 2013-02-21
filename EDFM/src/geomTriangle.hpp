/*!
 *	@file geomTriangle.hpp
 *	@brief Class for Triangle in 3D space.
 *
 *	@author Luca Turconi <lturconi@gmail.com>
 *  @date 22-09-2012
 *
 */

#ifndef GEOMTRIANGLE_HPP_
#define GEOMTRIANGLE_HPP_ 

#include<iostream>
#include<string>

#include "TypeDefinition.hpp"
#include "geomPoint3D.hpp"
#include "geomSegment.hpp"


namespace Geometry
{

	/*!
		@class Triangle
		
		@author Luca Turconi <lturconi@gmail.com>
		
    	This class implements the concept of Triangle.
    	It stores the three vertexes: M_pA, M_pB, M_pC.
    	
    	The class provides methods to compute the basic geometrical properties,
		such as area, perimeter and normal vector.
		
		It is also possible to test if a given point:
		 - lies in the plane generated by the triangle
		 - is contained in the triangle
		
		It is also possible to test for the existence of intersections
		with a segment. A method to compute the intersection point is also available.
		
		The method showMe allows to print the class informations.
    	It is also possible to print the attributes of the triangle by using << operator.
    	
    	For 3D visualization use the exportVtk method.
		
    */
class Triangle
{
public:
	//! @name Constructor & Destructor
	//@{
		
	//! Empty constructor
	Triangle();
	
	//! Constructor, getting the vertexes
	/*!
	 * @param a The first point
	 * @param b The second point
	 * @param c The third point
	 */
	Triangle(const Point3D & a, const Point3D & b, const Point3D & c);
	
	//! Copy constructor
	/*!
	 * @param t The segment copied in the new object
	 */
	Triangle(const Triangle & t);
	
	//! Destructor
	virtual ~Triangle();
	
	//@}
	
	//! @name Get Methods
	//@{
		
	//! Get vertex A
	/*!
	 * @return The vertex A
	 */
	inline Point3D A() const { return M_pA; }
	
	//! Get vertex B
	/*!
	 * @return The vertex B
	 */
	inline Point3D B() const { return M_pB; }
	
	//! Get vertex C
	/*!
	 * @return The vertex C
	 */
	inline Point3D C() const { return M_pC; }
	
	//! Get triangle maximum dimension
	/*!
	 * @return The M_Lmax value
	 */
	inline Real Lmax() const { return M_Lmax; }
	
	//@}
	
	//! @name Set Methods
	//@{
		
	//! Set vertex A
	/*!
	 * @param p The new value for the vertex A
	 */
	inline void setA(const Point3D & p) { M_pA=p; }
	
	//! Set vertex B
	/*!
	 * @param p The new value for the vertex B
	 */
	inline void setB(const Point3D & p) { M_pB=p; }
	
	//! Set vertex C
	/*!
	 * @param p The new value for the vertex C
	 */
	inline void setC(const Point3D & p) { M_pC=p; }
	
	//! Set Lmax
	/*!
	 * Set the maximum length of the triangle
	 */
	void setLmax();
	
	//@}
	
	//! @name Methods
	//@{
	
	//! Perimeter
	/*!
	 * It computes the perimeter of the triangle.
	 * @return The perimeter
	 */
	inline Real perimeter() const
		{ return (M_pA-M_pB).norm() + (M_pA-M_pC).norm() + (M_pC-M_pB).norm(); }
	
	//! Area
	/*!
	 * It computes the area of the triangle.
	 * @return The area
	 */
	inline Real area() const
		{ return 1./2.*((M_pB-M_pA).cross(M_pC-M_pA)).norm(); }
	
	//! Normal vector
	/*!
	 * It computes the normalized vector normal to the triangle.
	 * @return The normalized normal vector
	 */
	inline Vector3D normal() const
		{ return (M_pA-M_pB).cross(M_pA-M_pC) / ((M_pA-M_pB).cross(M_pA-M_pC).norm()); }
	
	//! Test coplanarity with a point
	/*!
	 * It tests if a given point lies on the plane generated by the triangle.
	 * @param p The point to be tested
	 * @return TRUE if the point is coplanar with the triangle, FALSE otherwise
	 */
	inline bool coplanarWithPoint(const Point3D & p) const
		{ return ( std::fabs(this->normal().dot(p-M_pA))
					<= eps*M_Lmax ); }
	
	//! Test if triangle contains a point
	/*!
	 * It tests if a given point belongs to the triangle.
	 * @param p The point to be tested
	 * @return TRUE if the point belongs to the triangle, FALSE otherwise
	 */
	bool containPoint(const Point3D & p) const;
	
	//! Test intersection with a segment
	/*!
	 * It tests if the triangle intersect a given segment.
	 * In the case of coplanar segment, the method will return FALSE.
	 * (A single intersection does not exist!)
	 * @param s The segment to be tested
	 * @return TRUE if the segment intersect the triangle
	 * 		FALSE if the intersection does not exist
	 */	
	bool isIntersectedBy(const Segment & s) const;
	
	//! Compute intersection with a segment
	/*!
	 * It compute the intersection between the triangle and a given segment.
	 * If the intersection does not exist, this method will return the point (NaN,NaN,NaN).
	 * Equally, in the case of coplanar segment, the method will return (NaN,NaN,NaN).
	 * (A single intersection does not exist!)
	 * @param s The segment
	 * @return The intersection point or (NaN,NaN,NaN) if the intersection does not exist.
	 */	
	Point3D intersectionWith(const Segment & s) const;

	std::vector<Point3D> getGaussNodes();
	
	std::vector<Real> getGaussWeights();

	//! Export in vtk format
	/*!
	 * It generates the vtk file, for the 3D visualization of the triangle.
	 * (Use paraview to open the vtk file)
	 * @param filename The name of the vtk file created by this method
	 * @return TRUE -> operation ended correctly
				FALSE -> an error occurred
	 */	
	bool exportVtk(const std::string & filename) const;
	
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
	
	Real M_Lmax;
};

//! @name External Operators
//@{
	
//! The insertion operator
/*!
 * @param ostr The stream object on which the action is performed.
 * 			This is the first parameter of the global functions,
 * 			and represents the object to the left of the operator,
 * 			i.e. the object on which the extraction operation is performed.
 * @param t The triangle inserted on the stream.
 * @return The stream object on which the action is performed (ostr)
 */
std::ostream& operator<<(std::ostream & ostr, const Triangle & t);
	
//@}
	
	
} // namespace Geometry

#endif /* GEOMTRIANGLE_HPP_ */
