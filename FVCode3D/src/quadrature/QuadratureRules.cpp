/*!
 *	@file QuadratureRules.cpp
 *	@brief Classes for quadrature rules (definitions).
 */

#include "quadrature/QuadratureRules.hpp"

namespace Darcy
{

std::unique_ptr<QuadratureRule> CentroidQuadrature::clone() const
{
	return std::unique_ptr<QuadratureRule>(new CentroidQuadrature(*this));
}

Real CentroidQuadrature::apply(const std::vector<Generic_Point> & pointVector, const Real volume, const std::function<Real(Generic_Point)> & Integrand) const
{
	Real q_integral;
	Generic_Point sum;
	UInt N = pointVector.size();

	for (UInt j = 0; j < N; ++j)
		sum = sum + (pointVector[j])/N;

	q_integral = Integrand(sum);

	return q_integral * volume;
}

} // namespace Darcy
