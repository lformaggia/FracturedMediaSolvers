/*!
 * @file preconditioner.hpp
 * @brief These classes implement simple classes that perform the "inversion" of the preconditioner (definition version)
*/

#include <FVCode3D/preconditioner/preconditioner.hpp>
#include <FVCode3D/core/BasicType.hpp>
#include <Eigen/LU>
#include <unsupported/Eigen/SparseExtra>

namespace FVCode3D
{

Vector diagonal_preconditioner::solve(const Vector & r)
{
	Vector z(A.rows());                     // The preconditioned residual
	for(UInt i = 0; i<A.rows(); i++)        // Solving Pz=r
		z[i] = r[i]/A.coeff(i,i);
	return z;                               // Move semantic will move the Eigen vector
}

}
