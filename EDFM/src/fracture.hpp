 /*!
 *	@file geomFault.hpp
 *	@brief Class for Fault.
 *
 *	@author Luca Turconi <lturconi@gmail.com>
 *  @date 22-09-2012
 *
 */

#ifndef FRACTURE_HPP_
#define FRACTURE_HPP_

#include "TypeDefinition.hpp"
#include "geomPoint3D.hpp"
#include "geomBilinearSurface.hpp"
#include "geomCPgridElements.hpp"
#include "geomCPgrid.hpp"
#include "interCellIntersections.hpp"
#include "interGridIntersections.hpp"
#include "interGridEdgeMap.hpp"
#include "interGridIntersectionMap.hpp"



namespace Geometry
{      class CProp;
	/*!
		@class Fault
		
		@author Luca Turconi <lturconi@gmail.com>
		
    	This class implements the concept of Fault.
    	
    	A fault is essentially a bilinear surface.
    	It provides additional methods to compute intersections with Corner Point
    	cells and Corner Point grids.
    	
    */
class Fracture: public BilinearSurface {
public:
	struct IntFrac {
	  Segment SMax;
	  Point3D Normale;
	  std::vector<gmm::size_type> M_i;
	  std::vector<gmm::size_type> M_j;
	  std::vector<gmm::size_type> M_k;
	  std::vector<gmm::size_type> Which1;
	  std::vector<gmm::size_type> Which2;
	  std::vector<Point3D> puntiInt;
	  std::vector<Real> M_ax;
	  std::vector<Real> M_ay;
	  std::vector<Real> M_az;
	  std::vector<Real> M_dmedio;
 	 } ; 
	//! @name Constructor & Destructor
	//@{
		
	//! Empty constructor
	Fracture();
	
	//! Constructor, getting the extremal points
	/*!
	 * Attention: when creating the fault, a and c are considered as opposite points.
	 * @param a The first point
	 * @param b The second point
	 * @param c The third point
	 * @param d The fourth point
	 */
	Fracture(const Point3D & a, const Point3D & b, const Point3D & c, const Point3D & d);

	//! Copy constructor
	/*!
	 * @param b The bilinear surface copied in the new object
	 */

	Fracture(const BilinearSurface & b);

	inline std::vector<Segment> getSx() const {return M_Sx;}

	inline std::vector<Segment> getSy() const {return M_Sy;}

	inline std::vector<Segment> getSz() const {return M_Sz;}

	inline std::vector<gmm::size_type> M_ii() const {return M_ipos;}

	inline std::vector<gmm::size_type> M_jj() const {return M_jpos;}

	inline std::vector<gmm::size_type> M_kk() const {return M_kpos;}

	void setProperties(Real, Real, Real);

	void setGeoProp(CProp& );
	
	void setNormaltoEdges();

	void setHalfLength();

	void setHalfTransmFF();

	void setTransmFF();

	void setTransmFM();
		
	inline Real perm() const {return M_perm;}

	inline Real compr() const {return M_compr;}

	inline Real aperture() const {return M_aperture;}
	
	//! Destructor
	virtual ~Fracture();

	inline void setIsInt(gmm::size_type i) {M_isintby.push_back(i);}

	inline void clearIsInt() {M_isintby.resize(0);}

	void setInt(gmm::size_type, Fracture &, std::vector<Point3D>, bool );

	inline std::vector<gmm::size_type> getIsInt() const {return M_isintby;}

	void computeNe();

	bool exportFracture(std::ofstream & , gmm::size_type);

	inline void setMetric(bool isMetric) {M_isMetric=isMetric;}

	inline std::vector<IntFrac> inter() {return M_inter;}

	Real CDARCY;
	
private:
	CPgrid* M_gridpointer; 
	bool M_isMetric;
	std::vector<gmm::size_type> M_isintby;	
std::vector<IntFrac> M_inter;
	Real M_perm;
	Real M_aperture;
	Real M_compr;
	gmm::size_type M_Ne;
	std::vector<gmm::size_type> M_ipos;
	std::vector<gmm::size_type> M_jpos;
	std::vector<gmm::size_type> M_kpos;
	std::vector<Real> M_areas;
	std::vector<Real> M_Dmedio;
	std::vector<Point3D> M_CG;
	std::vector<Segment> M_Sx;
	std::vector<Segment> M_Sy;
	std::vector<Segment> M_Sz;
	std::vector<Point3D> M_Lx;
	std::vector<Point3D> M_Ly;
	std::vector<Point3D> M_Lz;
	std::vector<Point3D> M_Nx;
	std::vector<Point3D> M_Ny;
	std::vector<Point3D> M_Nz;
	std::vector<Real> M_THX;
	std::vector<Real> M_THY;
	std::vector<Real> M_THZ;
	std::vector<Real> M_TX;
	std::vector<Real> M_TY;
	std::vector<Real> M_TZ;
	std::vector<Real> M_TM;

};

Segment maxSegment(std::vector<Point3D> &); 
//Real DmedioFF(std::vector<Point3D> &, Segment &,Fracture* );
	
} // namespace Geometry

#endif /* GEOMFAULT_HPP_ */
