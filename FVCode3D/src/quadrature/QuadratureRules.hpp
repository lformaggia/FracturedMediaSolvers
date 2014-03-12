/*!
 *	@file QuadratureRules.hpp
 *	@brief Classes for quadrature rules.
 */

#ifndef __DARCYQUADRATURERULE_HPP__
#define __DARCYQUADRATURERULE_HPP__

#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>

#include "core/TypeDefinition.hpp"

class Rigid_Mesh;

namespace Darcy{

//! Class that implements a quadrature rule
/*!
	@class QuadratureRule
	This is a base abstract class that implements a quadrature rule.
*/
class QuadratureRule
{
public:

	typedef Geometry::Point3D Generic_Point;
	typedef Geometry::Rigid_Mesh::Cell Cell;
	typedef Geometry::Rigid_Mesh::Facet Facet;
	typedef std::unique_ptr<QuadratureRule> QuadratureRuleHandler;

	//! Clone method for the class
	/*!
	 * It is used to pass the QuadratureRule to a class Quadrature
	 */
	virtual std::unique_ptr<QuadratureRule> clone() const = 0;

	//! Method apply for cells
	/*!
	@param cell the cell to integrate over
	@param Integrand The function we want to integrate
	@return the approximation of the integral of Integrand on the cell
	 */
	virtual Real apply(const Cell & cell, const std::function<Real(Generic_Point)> & Integrand) const = 0;

	//! Method apply for facets
	/*!
	@param facet the facet to integrate over
	@param volume volume of the facet
	@param Integrand The function we want to integrate
	@return the approximation of the integral of Integrand on the facet
	 */
	virtual Real apply(const Facet & facet, const Real volume, const std::function<Real(Generic_Point)> & Integrand) const = 0;

	//! Destructor
	virtual ~QuadratureRule() {};
};

//! MidPoint Quadrature rule for a generic polyhedron/polygon
/*!
	@class CentroidQuadrature
	This class is a quadrature rule for an arbitrary polyhedron/polygon
*/
class CentroidQuadrature : public QuadratureRule
{
public:

	//! Constructor
	CentroidQuadrature() = default;

	//! Clone method for the class
	/*!
	 * It is used to pass the QuadratureRule to a class Quadrature
	 */
    std::unique_ptr<QuadratureRule> clone() const;

	//! Method apply for cells
	/*!
	@param cell the cell to integrate over
	@param Integrand The function we want to integrate
	@return the approximation of the integral of Integrand on the cell
	 */
	virtual Real apply(const Cell & cell, const std::function<Real(Generic_Point)> & Integrand) const;

	//! Method apply for facets
	/*!
	@param facet the facet to integrate over
	@param volume volume of the facet
	@param Integrand The function we want to integrate
	@return the approximation of the integral of Integrand on the facet
	 */
	virtual Real apply(const Facet & facet, const Real volume, const std::function<Real(Generic_Point)> & Integrand) const;

	//! Destructor
    virtual ~CentroidQuadrature(){};

};

} // namespace Darcy

#endif
