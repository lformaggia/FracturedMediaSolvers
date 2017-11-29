/*!
 * @file preconditioner.hpp
 * @brief These classes implement simple classes that perform the "inversion" of the preconditioner (definition version)
*/

#include <FVCode3D/preconditioner/preconditioner.hpp>
#include <FVCode3D/assembler/global_operator.hpp>
#include <FVCode3D/assembler/local_operator.hpp>
#include <FVCode3D/core/BasicType.hpp>
#include <FVCode3D/core/Chrono.hpp>
#include <Eigen/LU>
#include <unsupported/Eigen/SparseExtra>

namespace FVCode3D
{

Vector diagonal_preconditioner::solve(const Vector & r) const
{
	Vector z = Vector::Zero(M.rows()+T.rows());      // The preconditioned residual
	for(UInt i = 0; i<M.rows(); i++)                 
	{
		if(M.coeff(i,i)!=0)
			z[i] = r[i]/M.coeff(i,i);
		else
			z[i] = r[i];
    }
	for(UInt i = 0; i<T.rows(); i++)                
	{
		if(T.coeff(i,i)!=0)
			z[M.rows()+i] = r[M.rows()+i]/T.coeff(i,i);          
		else
			z[M.rows()+i] = r[M.rows()+i];
    }    
	return z;                           
}

void lumpIP_builder::build(DiagMat & M_lump) const
{
	for (UInt k=0; k < M.outerSize(); ++k)
		for (SpMat::InnerIterator it(M,k); it; ++it)
			M_lump.diagonal()[it.row()] += it.value();
}

Vector BlockTriangular_preconditioner::solve(const Vector & r) const
{
	// First step: solve Inexact Schur Complement linear system
    Eigen::ConjugateGradient<SpMat> cg;
    cg.setMaxIterations(MaxIt);
    cg.setTolerance(tol);
    cg.compute(-ISC);
    Vector y2 = cg.solve(r.segment(Md_inv.rows(),B.rows()));
    // Second step: solve the diagonal linear system
    Vector z(Md_inv.rows()+B.rows());
	z.segment(0,Md_inv.rows()) = Md_inv*(r.segment(0,Md_inv.rows())+B.transpose()*y2);
    z.segment(Md_inv.rows(),B.rows()) = -y2;
    return z;
}

Vector ILU_preconditioner::solve(const Vector & r) const
{
	// First step: solve the 1st diagonal linear system
	Vector y1 = Md_inv*r.segment(0,Md_inv.rows());
	// Second step: solve the SC linear system
    Eigen::ConjugateGradient<SpMat> cg;
    cg.setMaxIterations(MaxIt);
    cg.setTolerance(tol);
    cg.compute(-ISC);
    Vector y2 = cg.solve(B*y1-r.segment(Md_inv.rows(),B.rows()));
    // Third step: solve the 2nd diagonal linear system
    Vector z(Md_inv.rows()+B.rows());
	z.segment(0,Md_inv.rows()) = y1-Md_inv*B.transpose()*y2;
    z.segment(Md_inv.rows(),B.rows()) = y2;
    return z;
}

}
