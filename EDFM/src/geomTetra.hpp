/*!
 *	@file geomTriangle.hpp
 *	@brief Class for Triangle in 3D space.
 *
 *	@author Luca Turconi <lturconi@gmail.com>
 *  @date 22-09-2012
 *
 */

#ifndef GEOMTETRA_HPP_
#define GEOMTETRA_HPP_ 

#include<iostream>
#include<string>

#include "TypeDefinition.hpp"
#include "geomPoint3D.hpp"
#include "geomSegment.hpp"
#include "geomTriangle.hpp"
#include "gmm/gmm.h"

namespace Geometry
{

	/*!
		@class Tetra
		
		@author Anna Scotti
	
		
    */
class Tetra
{
public:
	//! @name Constructor & Destructor
	//@{
		
	//! Empty constructor
	Tetra();
	
	Tetra(const Point3D & a, const Point3D & b, const Point3D & c, const Point3D & d);

	//! Copy constructor
	/*!
	 * @param t The segment copied in the new object
	 */
	Tetra(const Tetra & t);
	
	//! Destructor
	virtual ~Tetra();
	
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
	
	//! Get vertex D
	/*!
	 * @return The vertex D
	 */
	inline Point3D D() const { return M_pD; }

	Point3D getPoint(const gmm::size_type &);

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
	
	//! Set vertex D
	/*!
	 * @param p The new value for the vertex D
	 */
	inline void setD(const Point3D & p) { M_pD=p; }
	
	Triangle getFace(const Point3D &, const Point3D & , const Point3D &);
	
	Triangle getFace(const gmm::size_type &, const gmm::size_type &,const gmm::size_type &);
		
	//! Volume
	/*!
	 * It computes the area of the triangle.
	 * @return The area
	 */
	inline Real volume() const
		{ return 1./6.*fabs(((M_pB-M_pA).cross(M_pC-M_pA)).dot(M_pD-M_pA)); }
	
	std::vector<Point3D> getGaussNodes();

	std::vector<Real> getGaussWeights();
	
private:
	Point3D M_pA;
	Point3D M_pB;
	Point3D M_pC;
	Point3D M_pD;
};

	
	
} // namespace Geometry

#endif /* GEOMTRIANGLE_HPP_ */
