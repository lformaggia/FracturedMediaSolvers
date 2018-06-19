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
#include <chrono>
#include <cmath>

namespace FVCode3D
{

Vector diagonal_preconditioner::solve(const Vector & r) const
{
	auto & M = *Mptr;
	auto & T = *Tptr;
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

Vector BlockTriangular_preconditioner::solve(const Vector & r) const
{
	auto & B = *Bptr;
	// First step: solve Inexact Schur Complement linear system
	Vector y2 = chol.solve(r.tail(B.rows()));
    // Second step: solve the diagonal linear system
    Vector z(Md_inv.rows()+B.rows());
	z.head(Md_inv.rows()) = Md_inv*(r.head(Md_inv.rows())+B.transpose()*y2);
    z.tail(B.rows()) = -y2;
    return z;
}

Vector ILU_preconditioner::solve(const Vector & r) const
{
	auto & B = *Bptr;
	// First step: solve the 1st diagonal linear system
	Vector y1 = Md_inv*r.head(Md_inv.rows());
	// Second step: solve the ISC linear system
	Vector y2 = chol.solve(B*y1-r.tail(B.rows()));
    // Third step: solve the 2nd diagonal linear system
    Vector z(Md_inv.rows()+B.rows());
	z.head(Md_inv.rows()) = y1-Md_inv*B.transpose()*y2;
    z.tail(B.rows()) = y2;
    return z;
}

UInt constexpr HSS_preconditioner::MaxIt_default;
Real constexpr HSS_preconditioner::tol_default;
Real constexpr HSS_preconditioner::alpha_default;

void HSS_preconditioner::set(const SaddlePointMat & SP)
{
	alpha = 1.e-4;
	
	auto & M = SP.getM();
	auto & B = SP.getB();
	auto & T = SP.getT();
	Bptr = & SP.getB();
	Halpha = M;
	for(UInt i=0; i<M.rows(); i++)
		Halpha.coeffRef(i,i) += alpha;
		
	UInt c  = 0;   //To count number of fracture facets
	UInt cc = 0;
	for (UInt k=0; k<T.outerSize(); ++k)
	{
		for (SpMat::InnerIterator it(T,k); it; ++it)
		{
			if(it.value() != 0)
			{
				cc++;
				if(cc == 1)
					c++;
			}
		}
		cc = 0;
	}
	Nfrac = c;
	Ncell = T.rows()-c;
	SpMat Talpha = -T.bottomRightCorner(Nfrac,Nfrac);
	for(UInt i=0; i<Talpha.rows(); i++)
		Talpha.coeffRef(i,i) += alpha;		
		
	SpMat BBtalpha = B*B.transpose();
	for(UInt i=0; i<BBtalpha.rows(); i++)
		BBtalpha.coeffRef(i,i) += alpha*alpha;	
		
	cg.setMaxIterations(MaxIt);
	cg.setTolerance(tol);
	cg.compute(Halpha);
	cholT.compute(Talpha);
	cholBBt.compute(BBtalpha);
}

Vector HSS_preconditioner::solve(const Vector & r) const
{
	auto & B = *Bptr;
	// First step: solve the H linear system
	Vector omega1 = cg.solve(r.head(Halpha.rows()));
	// Second step: solve the T linear system
	Vector omega2(Ncell+Nfrac);
	for(UInt i = 0; i<Ncell; i++)                 
	{
		omega2[i] = r[Halpha.rows()+i]/alpha;
	}
	omega2.tail(Nfrac) = cholT.solve(r.tail(Nfrac));
    // Third step: solve the BBt linear system
    Vector z(Halpha.rows()+Ncell+Nfrac);
	z.tail(Ncell+Nfrac) = cholBBt.solve(B*omega1+alpha*omega2);
    z.head(Halpha.rows()) = (omega1-B.transpose()*z.tail(Ncell+Nfrac))/alpha;

    return z;
}

}
