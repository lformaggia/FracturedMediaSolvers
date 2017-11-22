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

Vector diagonal_preconditioner::solve(const Vector & r) const
{
	Vector z = Vector::Zero(A.cols());      // The preconditioned residual
	for(UInt i = 0; i<A.rows(); i++)        // Solving Pz=r
	{
		if(A.coeff(i,i)!=0)
			z[i] = r[i]/A.coeff(i,i);
		else
			z[i] = r[i];
    }
	return z;                               // Move semantic will move the Eigen vector
}

void lumpIP_builder::build(DiagMat & M_lump) const
{
	for (UInt k=0; k < M.outerSize(); ++k)
		for (SpMat::InnerIterator it(M,k); it; ++it)
			M_lump.diagonal()[it.row()] += it.value();
}

Vector BlockTriangular_preconditioner::solve(const Vector & r) const
{
	// First step: solve Inexact Schur Complement system
	Vector y2(B.rows());
	Eigen::SimplicialCholesky<SpMat, Eigen::Upper> chol(-ISC);
    y2 = chol.solve(r.segment(IM.rows(),B.rows()));
    // Second step: solve the Lumped system
    Vector z(IM.rows()+B.rows());
    z.segment(0,IM.rows()) = IM.inverse()*(r.segment(0,IM.rows())+B.transpose()*y2);
    z.segment(IM.rows(),B.rows()) = -y2;
    return z;
}


}
