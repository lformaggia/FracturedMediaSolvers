/*!
 *  @file QuadratureRules.hpp
 *  @brief Classes for quadrature rules.
 */

#ifndef __DARCYQUADRATURERULE_HPP__
#define __DARCYQUADRATURERULE_HPP__

#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>

#include <FVCode3D/core/TypeDefinition.hpp>

namespace FVCode3D
{

class Rigid_Mesh;

//! Class that implements a quadrature rule
/*!
 * @class QuadratureRule
 * This is a base abstract class that implements a quadrature rule.
 */
class QuadratureRule
{
public:

    typedef Rigid_Mesh::Cell Cell;
    typedef Rigid_Mesh::Facet Facet;
    typedef std::unique_ptr<QuadratureRule> QuadratureRuleHandler;

    //! Constructor
    QuadratureRule() = default;

    //! Clone method for the class
    /*!
     * It is used to pass the QuadratureRule to a class Quadrature
     */
    virtual std::unique_ptr<QuadratureRule> clone() const = 0;

    //! Method apply for cells
    /*!
     * @param cell the cell to integrate over
     * @param Integrand The function we want to integrate
     * @return the approximation of the integral of Integrand on the cell
     */
    virtual Real apply(const Cell & cell, const std::function<Real(Point3D)> & integrand) const = 0;

    //! Method apply for facets
    /*!
     * @param facet the facet to integrate over
     * @param volume volume of the facet
     * @param Integrand The function we want to integrate
     * @return the approximation of the integral of Integrand on the facet
     */
    virtual Real apply(const Facet & facet, const std::function<Real(Point3D)> & integrand) const = 0;

    //! Destructor
    virtual ~QuadratureRule() = default;
};

//! MidPoint Quadrature rule for a generic polyhedron/polygon
/*!
 * @class CentroidQuadrature
 * This class is a quadrature rule for an arbitrary polyhedron/polygon
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
     * @param cell the cell to integrate over
     * @param Integrand The function we want to integrate
     * @return the approximation of the integral of Integrand on the cell
     */
    virtual Real apply(const Cell & cell, const std::function<Real(Point3D)> & integrand) const;

    //! Method apply for facets
    /*!
     * @param facet the facet to integrate over
     * @param volume volume of the facet
     * @param Integrand The function we want to integrate
     * @return the approximation of the integral of Integrand on the facet
     */
    virtual Real apply(const Facet & facet, const std::function<Real(Point3D)> & integrand) const;

    //! Destructor
    virtual ~CentroidQuadrature() = default;
};

//! Tetrahedralization quadrature rule for a generic polyhedron
/*!
 * @class Tetrahedralization
 * The polyhedron is decomposed into tetrahedras and over each tetrahedra a midpoint quadrature rule is used.
 * The facet version of this method builds up a triangulation of the polygon to compute the integral.
 */
class Tetrahedralization : public QuadratureRule
{
public:

    //! Constructor
    Tetrahedralization() = default;

    //! Clone method for the class
    /*!
     * It is used to pass the QuadratureRule to a class Quadrature
     */
    std::unique_ptr<QuadratureRule> clone() const;

//    void TetFunc(const std::vector<Point3D> & tetrahedronNodes, std::vector<std::vector<Point3D>> & V) const;
    
//    void TetFunc2(const std::vector<Point3D> & tetrahedronNodes, const UInt n, Real & result,
//		const std::function<Real(Point3D)> & integrand) const;

    //! Method apply for cells
    /*!
     * @param cell the cell to integrate over
     * @param Integrand The function we want to integrate
     * @return the approximation of the integral of Integrand on the cell
     */
    virtual Real apply(const Cell & cell, const std::function<Real(Point3D)> & integrand) const;

    //! Method apply for facets
    /*!
     * @param facet the facet to integrate over
     * @param Integrand The function we want to integrate
     * @return the approximation of the integral of Integrand on the facet
     */
    virtual Real apply(const Facet & facet, const std::function<Real(Point3D)> & integrand) const;

    //! Destructor
    virtual ~Tetrahedralization() = default;
};

} // namespace FVCode3D
#endif
