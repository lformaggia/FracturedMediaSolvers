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
	  gmm::size_type altra;
	  Segment SMax;
	  Point3D Normale;
	  std::vector<gmm::size_type> M_i;
	  std::vector<gmm::size_type> M_j;
	  std::vector<gmm::size_type> M_k;
	  std::vector<gmm::size_type> Which1;
	  std::vector<gmm::size_type> Which2;
	  std::vector<gmm::size_type> Which1sorted;
	  std::vector<gmm::size_type> Which2sorted;
	  std::vector<Point3D> puntiInt;
	  std::vector<Real> M_ax;
	  std::vector<Real> M_ay;
	  std::vector<Real> M_az;
	  std::vector<Real> M_dmedio;
	  std::vector<Real> M_Tf1f2h;
	  std::vector<Real> M_Tf1f2;

 	 } ; 
	struct CellaF{     
	std::vector<Point3D> puntiAree;
	std::vector<Point3D> raggi;
	std::vector<Real> theta;
	};
	struct Neigh_cell{
	gmm::size_type nnc;
	std::vector<gmm::size_type> vicini;
	std::vector<Real> TFF;
	};

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

	inline std::vector<Segment> getSx(gmm::size_type i) const {return (i==1)?M_S1x:M_S2x;}

	inline std::vector<Segment> getSy(gmm::size_type i) const {return (i==1)?M_S1y:M_S2y;}

	inline std::vector<Segment> getSz(gmm::size_type i) const {return (i==1)?M_S1z:M_S2z;}

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

	//bool exportFracture(std::ofstream & , gmm::size_type);
	bool exportFracture2(std::ofstream & , gmm::size_type);

	inline void setMetric(bool isMetric) {M_isMetric=isMetric; CDARCY=(M_isMetric==true)?0.008527:0.001127;}

	inline std::vector<IntFrac> inter() {return M_inter;}

	bool isAreaOK() ;

	bool exportVtk(const std::string & ) const;
	
	void setInterTransm(gmm::size_type quale, gmm::size_type dove, Real quanto) {M_inter[quale].M_Tf1f2[dove]=quanto;}

	Real CDARCY;

	void sortCG_Y(); 	
	void sortCG_X(); 	
	std::vector<gmm::size_type> getCGsort() {return M_sortCG;}
	
private:
	CPgrid* M_gridpointer; 
	bool M_isMetric;
	std::vector<gmm::size_type> M_isintby;	
	std::vector<IntFrac> M_inter;
	std::vector<CellaF> M_celle;  
	Real M_perm;
	Real M_aperture;
	Real M_compr;
	gmm::size_type M_Ne;
	std::vector<gmm::size_type> M_ipos;
	std::vector<gmm::size_type> M_jpos;
	std::vector<gmm::size_type> M_kpos;
	std::vector<gmm::size_type> M_sortCG;
	std::vector<Neigh_cell> M_vicine;
	std::vector<Neigh_cell> M_vicineH;
	std::vector<Neigh_cell> M_vicineV;
	std::vector<Real> M_areas;
	std::vector<Real> M_Dmedio;
	std::vector<Point3D> M_CG;
	std::vector<Segment> M_S1x;
	std::vector<Segment> M_S1y;
	std::vector<Segment> M_S1z;
	std::vector<Segment> M_S2x;
	std::vector<Segment> M_S2y;
	std::vector<Segment> M_S2z;

	std::vector<Point3D> M_L1x;
	std::vector<Point3D> M_L1y;
	std::vector<Point3D> M_L1z;
	std::vector<Point3D> M_L2x;
	std::vector<Point3D> M_L2y;
	std::vector<Point3D> M_L2z;
	std::vector<Point3D> M_N1x;
	std::vector<Point3D> M_N1y;
	std::vector<Point3D> M_N1z;
	std::vector<Point3D> M_N2x;
	std::vector<Point3D> M_N2y;
	std::vector<Point3D> M_N2z;

	std::vector<Real> M_TH1X;
	std::vector<Real> M_TH1Y;
	std::vector<Real> M_TH1Z;
	std::vector<Real> M_TH2X;
	std::vector<Real> M_TH2Y;
	std::vector<Real> M_TH2Z;
	std::vector<Real> M_TX;
	std::vector<Real> M_TY;
	std::vector<Real> M_TZ;
	std::vector<Real> M_TM;

};

Segment maxSegment(std::vector<Point3D> &); 
//Real DmedioFF(std::vector<Point3D> &, Segment &,Fracture* );
	
} // namespace Geometry

#endif /* GEOMFAULT_HPP_ */
