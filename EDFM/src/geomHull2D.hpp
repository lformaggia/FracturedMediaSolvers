/*!
 *	@file geomTriangle.hpp
 *	@brief Class for Triangle in 3D space.
 *
 *	@author Luca Turconi <lturconi@gmail.com>
 *  @date 22-09-2012
 *
 */

#ifndef GEOMHULL2D_HPP_
#define GEOMHULL2D_HPP_ 

#include<iostream>
#include<string>

#include "TypeDefinition.hpp"
#include "geomPoint3D.hpp"
#include "geomTriangle.hpp"
#include "geomTetra.hpp"
#include "geomFault.hpp"
#include "gmm/gmm.h"

  extern "C" {
#ifdef _MSC_VER
# include <libqhull/qhull_a.h>
#else
# include <qhull/qhull.h>
//# include <qhull/mem.h>
# include <qhull/qset.h>
# include <qhull/geom.h>
# include <qhull/merge.h>
# include <qhull/poly.h>
# include <qhull/io.h>
# include <qhull/stat.h>
#endif
}
namespace Geometry
{

	/*!
		@class Hull
		
		@author Anna Scotti
		
    	*/
class Hull2D
{
public:

	//! @name Constructor & Destructor
	//@{
		
	//! Empty constructor
	Hull2D();
	
	//! Constructor, getting the simplexes matrix

	Hull2D(std::vector<Point3D>, Fault*);

	//! Destructor
	virtual ~Hull2D();
	
	//@}
	
	//! @name Get Methods
	//@{
		
	//! Get vertex A
	/*!
	 * @return The vertex A
	 */
	inline gmm::size_type getNtriangle() const {return M_Ntriangle;}
	inline Point3D getPoint(gmm::size_type & i) const { return M_points[i]; }
	inline gmm::dense_matrix<gmm::size_type> getSimplexes() const { return M_simplexes; }
	
	bool call_qhull2D(gmm::size_type & , gmm::dense_matrix<gmm::size_type> &, std::vector<coordT> &);

	inline std::vector<Triangle> getTriangle() {return M_triangle;}

	inline Triangle getTriangle(gmm::size_type i) {return M_triangle[i];}

	Real getArea();

	std::vector<gmm::size_type> getPointsSimplex(gmm::size_type );

private:
	gmm::dense_matrix<gmm::size_type> M_simplexes;
	std::vector<Point3D> M_points;
	std::vector<Triangle> M_triangle;
	gmm::size_type M_Ntriangle;
};

	
	
} // namespace Geometry

#endif /* GEOMTRIANGLE_HPP_ */
