 /*!
 *	@file geomFault.hpp
 *	@brief Class for Fault.
 *
 *	@author Luca Turconi <lturconi@gmail.com>
 *  @date 22-09-2012
 *
 */

#ifndef FRACTURES_HPP_
#define FRACTURES_HPP_

#include "TypeDefinition.hpp"
#include "geomPoint3D.hpp"
#include "geomBilinearSurface.hpp"
#include "fracture.hpp"
#include "geomCPgridElements.hpp"
#include "geomCPgrid.hpp"
#include "interCellIntersections.hpp"
#include "interGridIntersections.hpp"
#include "interGridEdgeMap.hpp"
#include "interGridIntersectionMap.hpp"
#include "gmm/gmm.h"

namespace Geometry
{

class Fractures {
public:
	//! @name Constructor & Destructor
	//@{
		
	//! Empty constructor
	Fractures();
	
	//! Constructor
	Fractures(const std::string , Real);
	
	//! Destructor
	virtual ~Fractures();
	
	void computeIntersections(bool );
	inline gmm::size_type getNfractures() const {return M_nfractures;}
	//inline std::vector<Fracture> getfractures() const {return M_fractures;}
	//inline Fracture getfracture(gmm::size_type i) const {return M_fractures[i];}

	void setInterFTransm();

	gmm::size_type M_nfractures;
	std::vector<Fracture> M_fractures;
	bool M_isMetric;
	
};

void readPoints(std::string line, Real &, Real &, Real &);
void readProperties(std::string line, gmm::size_type &, Real &, Real &, Real &);
} // namespace Geometry

#endif /* GEOMFAULT_HPP_ */
