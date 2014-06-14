/*!
 *	@file QuadratureRules.cpp
 *	@brief Classes for quadrature rules (definitions).
 */

#include "mesh/Rigid_Mesh.hpp"
#include "quadrature/QuadratureRules.hpp"

namespace FVCode3D
{

std::unique_ptr<QuadratureRule> CentroidQuadrature::clone() const
{
	return std::unique_ptr<QuadratureRule>(new CentroidQuadrature(*this));
}

Real CentroidQuadrature::apply(const Cell & cell, const std::function<Real(Point3D)> & integrand) const
{
	return integrand(cell.getCentroid()) * cell.getVolume();
}

Real CentroidQuadrature::apply(const Facet & facet, const Real volume, const std::function<Real(Point3D)> & integrand) const
{
	return integrand(facet.getCentroid()) * volume;
}

} // namespace FVCode3D
