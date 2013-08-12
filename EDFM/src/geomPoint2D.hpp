/*!
 *	@file geomPoint3D.hpp
 *	@brief Struct for Point and Vector in 3D space.
 *
 *	@author Luca Turconi <lturconi@gmail.com>
 *  @date 22-09-2012
 *
 */

#ifndef GEOMPOINT2D_HPP_
#define GEOMPOINT2D_HPP_

#include<iostream>
#include<cmath>

#include "TypeDefinition.hpp"


namespace Geometry
{
	
	/*!
		@struct Point2D
		
		@author Luca Formaggia
    	
    	This struct implements the concept of 2D point.
    	It stores the cartesian coordinates (x,y).
    	
    	Each point is considered as the corresponding position vector.
    	According to this point of view the code implements:
    	
    	- basic vector operations, such as dot and cross products;
		
		- an overloading of the operators: +, -, *, /.
		
		The method showMe allows to print the struct informations.
		It is also possible to print the point by using << operator.
    	
    */
struct Point2D{
	
	//! @name Constructor & Destructor
	//@{
		
	//! Empty constructor
	Point2D();
	
	//! Constructor, getting the coordinates
	/*!
	 * @param a The x coordinate
	 * @param b The y coordinate
	 */
	Point2D(const Real & a, const Real & b);
	
	//! Copy constructor
	/*!
	 * @param p The point copied in the new object
	 */
	Point2D(const Point2D & p);
	
	//! Destructor
	~Point2D();
	
	//@}
	
	
	//! @name Methods
	//@{
	
	//! Dot products
	/*!
	 * It performs the dot product between <b>this</b> and p.
	 * @param p The second arguments of the dot product
	 * @return The dot product result (a real number)
	 */
	inline Real dot(const Point2D & p) const
		{ return (x*p.x + y*p.y); }
	
	//! The Norm
	/*!
	 * It computes the norm of the position vector.
	 * @return The position vector norm
	 */
	inline Real norm() const
		{ return std::sqrt(this->dot(*this)); }
	
  //! Scale the point.
  /*!
    Points are scaled as \f$ (x-origin)/factors\f$
    @parameter orig The origin.
    @parameter factors. A Real[2] array with the scaling factors.
   */
  void   scale(Point2D const & orig, Real factors[]);
  //! Scale back the points.
  /*!
    Points are scaled back as \f$ origin + x*factors\f$
    @param orig The origin.
    @param factors. A Real[2] array with the scaling factors.
   */
  void   scaleBack(Point2D const & orig, Real factors[]);
	//! Cross products
	/*!
	 * It performs the cross product between <b>this</b> and p.
	 * @param p The second arguments of the cross product
	 * @return The cross product result (a Real)
	 */
	Real cross(const Point2D & p) const;
	
	//! Display general information about the content of the struct
	/*!
	 * List of things displayed in the struct
	 * @param out Specify the output format (std::cout by default)
	 */
	void showMe(std::ostream  & out=std::cout) const;
	
	//@}
	
	//! @name Operators
	//@{
		
	//! Add a point to <b>this</b>
	/*!
	 * @param b The added point
	 * @return The operation result (a point)
	 */
	bool operator==( const Point2D & b)const;

	Point2D & operator+=( const Point2D & b);
	
	//! Subtract a point to <b>this</b>
	/*!
	 * @param b The subtracted point
	 * @return The operation result (a point)
	 */
	Point2D & operator-=( const Point2D & b);
	
	//! Multiply <b>this</b> by a scalar
	/*!
	 * @param a The scalar coefficient
	 * @return The operation result (a point)
	 */
	Point2D & operator*=( const Real & a);
	
	//! Divide <b>this</b> by a scalar
	/*!
	 * @param a The scalar coefficient
	 * @return The operation result (a point)
	 */
	Point2D & operator/=( const Real & a);

	//@}
	
    Real x;
    Real y;
};

//! @name External Operators
//@{

//! The sum operator for two points
/*!
 * @param p1 The first addend
 * @param p2 The second addend
 * @return The operation result (a point)
 */
Point2D const operator+(const Point2D & p1, const Point2D & p2);

//! The difference operator for two points
/*!
 * @param p1 The first addend
 * @param p2 The second addend
 * @return The operation result (a point)
 */
Point2D const operator-(const Point2D & p1, const Point2D & p2);

//! The scalar multiplication operator
/*!
 * @param a The scalar
 * @param p The point
 * @return The operation result (a point)
 */
Point2D const operator*(const Real & a, const Point2D & p);

//! The scalar multiplication operator (commutative property)
/*!
 * @param p The point
 * @param a The scalar
 * @return The operation result (a point)
 */
Point2D const operator*(const Point2D & p, const Real & a);

//! The scalar division operator (commutative property)
/*!
 * @param p The point
 * @param a The scalar
 * @return The operation result (a point)
 */
Point2D const operator/(const Point2D & p, const Real & a);

//! The insertion operator
/*!
 * @param ostr The stream object on which the action is performed.
 * 			This is the first parameter of the global functions,
 * 			and represents the object to the left of the operator,
 * 			i.e. the object on which the extraction operation is performed.
 * @param p The point inserted on the stream.
 * @return The stream object on which the action is performed (ostr)
 */
std::ostream& operator<<(std::ostream & ostr, const Point2D & p);

bool operator<(const Point2D & A, const Point2D & B);

//@}

	/*!
		@typedef Vector2D
    	2D vectors have the same structure and the same properties of 2D points.
    */
typedef Point2D Vector2D;


} // namespace Geometry

#endif /* GEOMPOINT2D_HPP_ */
