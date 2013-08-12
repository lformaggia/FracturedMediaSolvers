/*!
 *	@file geomPoint3D.hpp
 *	@brief Struct for Point and Vector in 3D space.
 *
 *	@author Luca Turconi <lturconi@gmail.com>
 *  @date 22-09-2012
 *
 */

#ifndef GEOMPOINT3D_HPP_
#define GEOMPOINT3D_HPP_

#include<iostream>
#include<cmath>

#include "TypeDefinition.hpp"
#include "tolerances.hpp"

namespace Geometry
{
	
	/*!
		@struct Point3D
		
		@author Luca Turconi <lturconi@gmail.com>
    	
    	This struct implements the concept of 3D point.
    	It stores the cartesian coordinates (x,y,z).
    	
    	Each point is considered as the corresponding position vector.
    	According to this point of view the code implements:
    	
    	- basic vector operations, such as dot and cross products;
		
		- an overloading of the operators: +, -, *, /.
		
		The method showMe allows to print the struct informations.
		It is also possible to print the point by using << operator.
    	
    */
  struct Point3D{
    
    //! @name Constructor & Destructor
    //@{
    
    //! Empty constructor
    Point3D();
    
    //! Constructor, getting the coordinates
    /*!
     * @param a The x coordinate
     * @param b The y coordinate
     * @param c The z coordinate
     */
    Point3D(const Real & a, const Real & b, const Real & c);
    
    //! Copy constructor
    /*!
     * @param p The point copied in the new object
     */
    Point3D(const Point3D & p);
    
    //! Destructor
    ~Point3D();
    
    //@}
    
    
    //! @name Methods
	//@{
	
	//! Dot products
	/*!
	 * It performs the dot product between <b>this</b> and p.
	 * @param p The second arguments of the dot product
	 * @return The dot product result (a real number)
	 */
	inline Real dot(const Point3D & p) const
		{ return (x*p.x + y*p.y + z*p.z); }
	
	//! The Norm
	/*!
	 * It computes the norm of the position vector.
	 * @return The position vector norm
	 */
	inline Real norm() const
		{ return std::sqrt(this->dot(*this)); }
	
	//! Cross products
	/*!
	 * It performs the cross product between <b>this</b> and p.
	 * @param p The second arguments of the cross product
	 * @return The cross product result (a point)
	 */
	Point3D cross(const Point3D & p) const;
	
	//! Display general information about the content of the struct
	/*!
	 * List of things displayed in the struct
	 * @param out Specify the output format (std::cout by default)
	 */
	void showMe(std::ostream  & out=std::cout) const;
	
	//@}
	
	//! @name Operators
	//@{

        //! Equality operator
        /*!
	  It compares two points using a tolerance.
	  @param b Point to be compared.
	  @return true if points are nearer than the tolerance
	  @note It does not compy with the semantic of the equality operator. It should be
	  changed asap by introducing a method instead.
	 */
	bool operator==( const Point3D & b) const;

	//! Add a point to <b>this</b>
	/*!
	 * @param b The added point
	 * @return The operation result (a point)
	 */
	Point3D & operator+=( const Point3D & b);
	
	//! Subtract a point to <b>this</b>
	/*!
	 * @param b The subtracted point
	 * @return The operation result (a point)
	 */
	Point3D & operator-=( const Point3D & b);
	
	//! Multiply <b>this</b> by a scalar
	/*!
	 * @param a The scalar coefficient
	 * @return The operation result (a point)
	 */
	Point3D & operator*=( const Real & a);
	
	//! Divide <b>this</b> by a scalar
	/*!
	 * @param a The scalar coefficient
	 * @return The operation result (a point)
	 */
	Point3D & operator/=( const Real & a);

	//@}
	
    Real x;
    Real y;
    Real z;
};

Point3D apply_shear(Point3D p, const Real m, std::string direction);
//! @name External Operators
//@{

//! The sum operator for two points
/*!
 * @param p1 The first addend
 * @param p2 The second addend
 * @return The operation result (a point)
 */
Point3D const operator+(const Point3D & p1, const Point3D & p2);

//! The difference operator for two points
/*!
 * @param p1 The first addend
 * @param p2 The second addend
 * @return The operation result (a point)
 */
Point3D const operator-(const Point3D & p1, const Point3D & p2);

//! The scalar multiplication operator
/*!
 * @param a The scalar
 * @param p The point
 * @return The operation result (a point)
 */
Point3D const operator*(const Real & a, const Point3D & p);

//! The scalar multiplication operator (commutative property)
/*!
 * @param p The point
 * @param a The scalar
 * @return The operation result (a point)
 */
Point3D const operator*(const Point3D & p, const Real & a);

//! The scalar division operator (commutative property)
/*!
 * @param p The point
 * @param a The scalar
 * @return The operation result (a point)
 */
Point3D const operator/(const Point3D & p, const Real & a);

//! The insertion operator
/*!
 * @param ostr The stream object on which the action is performed.
 * 			This is the first parameter of the global functions,
 * 			and represents the object to the left of the operator,
 * 			i.e. the object on which the extraction operation is performed.
 * @param p The point inserted on the stream.
 * @return The stream object on which the action is performed (ostr)
 */
std::ostream& operator<<(std::ostream & ostr, const Point3D & p);

//! Less-than operator
/*!
 * implements a lexicografic almost ordering
 * @param A the first point
 * @param B the second point
 * @return true if A is "smaller" than B according to a quasi-lexicografic ordering.
 * @note This operator does not comply with the good ordering rule. It must be changed asap.
 */

bool operator<(const Point3D & A, const Point3D & B);

//@}

	/*!
		@typedef Vector3D
    	3D vectors have the same structure and the same properties of 3D points.
    */
typedef Point3D Vector3D;


} // namespace Geometry

#endif /* GEOMPOINT3D_HPP_ */
