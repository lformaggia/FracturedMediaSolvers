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

namespace Darcy{

//! Class that implements a quadrature rule
/*!
	@class QuadratureRule
	This is a base abstract class that implements a quadrature rule.
*/
class QuadratureRule
{
public:

	typedef typename Geometry::Point3D Generic_Point;
	typedef typename std::unique_ptr< QuadratureRule > QuadratureRuleHandler;

	//! Clone method for the class
	/*!
	 * It is used to pass the QuadratureRule to a class Quadrature
	 */
	virtual std::unique_ptr<QuadratureRule> clone() const =0;

	//! Method apply
	/*!
	@param pointVector A vector of the vertexes of a cell
	@param volume The volume of the cell
	@param Integrand The function we want to integrate
	@return the approximation of the integral of Integrand on the cell
	 */
	virtual Real apply(const std::vector<Generic_Point> & pointVector, const Real volume, const std::function<Real(Generic_Point)> & Integrand) const = 0;

	//! destructor
	virtual ~QuadratureRule(){};
};

//! MidPoint Quadrature rule for a generic polyhedron
/*!
	@class CentroidQuadrature
	This class is a quadrature rule for an arbitrary polyhedron
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
    std::unique_ptr<QuadratureRule > clone() const;

	//! Method apply
	/*!
	@param pointVector A vector of the vertexes of a cell
	@param volume The volume of the cell
	@param Integrand The function we want to integrate
	@return the approximation of the integral of Integrand on the cell
	 */
	Real apply (const std::vector<Generic_Point> & pointVector, const Real volume, const std::function<Real(Generic_Point)> & Integrand) const;

	//! Destructor
    virtual ~CentroidQuadrature(){};

};

} // namespace Darcy

#endif
