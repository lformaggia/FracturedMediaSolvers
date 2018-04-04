/*!
 *  @file Quadrature.hpp
 *  @brief Class for integrating functions on a Rigid_Mesh.
 */

#ifndef __DARCYQUADRATURE_HPP__
#define __DARCYQUADRATURE_HPP__

#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>

#include <FVCode3D/core/TypeDefinition.hpp>
#include <FVCode3D/quadrature/QuadratureRules.hpp>

namespace FVCode3D
{

class Rigid_mesh;
class PropertiesMap;

//! Class that implements method to integrate function over the mesh
/*!
 * @class Quadrature
 * This class integrates solutions to PDE solved with the finite volume method,
 * computes the L2Norm of such a solution and integrates continuous functions.
 * By passing a QuadratureRule to the constructor it is possible to choose the quadrature rule you want to use.
 */
class Quadrature
{
    //! Typedef for QuadratureRule::QuadratureRuleHandler
    /*!
     * @typedef QR_Handler
     * This type definition permits to handle a QuadratureRule::QuadratureRuleHandler as a QR_Handler.
     */
    typedef QuadratureRule::QuadratureRuleHandler QR_Handler;

public:
    //! @name Constructor & Destructor
    //@{

    //! Constructor for Quadrature, given a Rigid_Mesh.
    /*!
     * @param rigid_mesh A Rigid_Mesh on which we want to integrate
     */
    Quadrature (const Rigid_Mesh & rigid_mesh);

    //! Constructor for Quadrature, given a Rigid_Mesh and a QuadratureRule
    /*!
     * @param rigid_mesh A Rigid_Mesh on which we want to integrate
     * @param quadrature A QuadratureRule we want to use in order to integrate a function
     */
    Quadrature (const Rigid_Mesh & rigid_mesh, const QuadratureRule & quadrature);

    //! Constructor for Quadrature, given a Rigid_Mesh, a QuadratureRule, and a QuadratureRule to integrate fractures
    /*!
     * @param rigid_mesh A Rigid_Mesh on which we want to integrate
     * @param quadrature A QuadratureRule we want to use in order to integrate a function
     * @param fracturequadrature A QuadratureRule we want to use in order to integrate a function on fractures
     */
    Quadrature (const Rigid_Mesh & rigid_mesh, const QuadratureRule & quadrature, const QuadratureRule & fracturequadrature);

    //! Destructor
    ~Quadrature() = default;

    //@}

    //! @name Methods
    //@{

    //! Integrate a function
    /*!
     * @param Integrand A std::function which returns a scalar
     * @return The integral of the considered Integrand function
     */
    Real integrate (const std::function<Real(Point3D)> & integrand) const;

    //! Integrate a function only on the porous matrix
    /*!
     * @param Integrand A std::function which returns a scalar
     * @return The integral of the considered Integrand function
     */
    Real integrateMatrix (const std::function<Real(Point3D)> & integrand) const;

    //! Integrate a function only on the fractures
    /*!
     * @param Integrand A std::function which returns a scalar
     * @return The integral of the considered Integrand function
     */
    Real integrateFractures (const std::function<Real(Point3D)> & integrand) const;


    //! Integrate discrete function
    /*!
     * @param Integrand A vector such that in the i-th component has the value of the function on the i-th cell
     * @return The integral of the considered Integrand function
     */
    Real integrate (const Vector & integrand) const throw();
    //! Integrate discrete function only on the porous matrix
    /*!
     * @param Integrand A vector such that in the i-th component has the value of the function on the i-th cell
     * @return The integral of the considered Integrand function
     */
    Real integrateMatrix (const Vector & integrand) const throw();
    //! Integrate discrete function only on the fracturs
    /*!
     * @param Integrand A vector such that in the i-th component has the value of the function on the i-th cell
     * @return The integral of the considered Integrand function
     */
    Real integrateFractures (const Vector & integrand) const throw();


    //! L2 Norm of a discrete function
    /*!
     * @param Integrand A vector such that in the i-th component has the value of the function on the i-th cell
     * @return The L2 norm of the considered Integrand function
     */
    Real L2Norm (const Vector & integrand) const;

    //! L2 Norm of a discrete function, only on the porous matrix
    /*!
     * @param Integrand A vector such that in the i-th component has the value of the function on the i-th cell
     * @return The L2 norm of the considered Integrand function
     */
    Real L2NormMatrix (const Vector & integrand) const;

    //! L2 Norm of a discrete function, only on the fractures
    /*!
     * @param Integrand A vector such that in the i-th component has the value of the function on the i-th cell
     * @return The L2 norm of the considered Integrand function
     */
    Real L2NormFractures (const Vector & integrand) const;


    //! Integrate a function and return the integral cell by cell
    /*!
     * @param Integrand A function which returns a scalar
     * @return A vector with in the i-th component the integral of the considered Integrand function on the i-th cell
     */
    Vector cellIntegrate (const std::function<Real(Point3D)> & func) const;

    //! Integrate a function and return the integral cell by cell, only in the porous matrix
    /*!
     * @param Integrand A function which returns a scalar
     * @return A vector with in the i-th component the integral of the considered Integrand function on the i-th cell
     * of the porous matrix
     * @note The vector contains all the problem entries, fractures included which are zero
     */
    Vector cellIntegrateMatrix (const std::function<Real(Point3D)> & func) const;

    //! Integrate a function and return the integral cell by cell, only in the fractures
    /*!
     * @param Integrand A function which returns a scalar
     * @return A vector with in the i-th component the integral of the considered Integrand function on the i-th cell,
     * of the fractures
     * @note The vector contains all the problem entries, the entries of the matrix are zero
     */
    Vector cellIntegrateFractures (const std::function<Real(Point3D)> & func) const;

    //@}

    //! @name Get Methods
    //@{

    //! Get Mesh size (const)
    /*!
     * @return The number of the cells in the Rigid_Mesh considered
     */
    UInt getMeshSize() const
        {return M_size;}
    //@}

protected:
    //! A reference to a Rigid_Mesh
    const Rigid_Mesh & M_mesh;
    //! A reference to a PropertiesMap
    const PropertiesMap & M_properties;
    //! The number of cells in M_mesh
    UInt M_size;
    //! A pointer to a QuadratureRule
    QR_Handler M_quadrature;
    //! A pointer to a QuadratureRule which is used to integrate on fractures
    QR_Handler M_fractureQuadrature;
};

} // namespace FVCode3D
#endif
