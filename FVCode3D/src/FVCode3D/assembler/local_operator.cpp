/*!
 * @file local_operator.cpp
 * @brief These classes implement local mimetic operators (definitions).
 */

#include <vector>
#include <cmath>
#include <FVCode3D/assembler/local_operator.hpp>
#include <FVCode3D/mesh/RigidMesh.hpp>
#include <FVCode3D/core/BasicType.hpp>
#include <Eigen/LU>
#include <unsupported/Eigen/SparseExtra>

namespace FVCode3D
{
	
constexpr Real local_InnerProduct::gamma;

void local_InnerProduct::assemble()      
{	
	Mat      Sp;         // I will take Np as the column base of Sp
	Mat33    Kp;         // permability tensor
	
	auto & K = pMesh.getPropertiesMap().getProperties(cellp.getZoneCode()).M_permeability; 
	Real Mob = pMesh.getPropertiesMap().getMobility();
	const std::vector<UInt> & cellFacetsId( cellp.getFacetsIds() );
	const UInt numCellFacets  = cellFacetsId.size();
    const Real cellMeasure    = cellp.getVolume();

    // Resize local matrices
    Np.resize(numCellFacets,Eigen::NoChange);
	Np.setZero();
    Rp.resize(numCellFacets,Eigen::NoChange);
	Rp.setZero();

	Kp(0,0) = Mob*K->operator()(0,0);
	Kp(0,1) = Mob*K->operator()(0,1);
	Kp(0,2) = Mob*K->operator()(0,2);
	Kp(1,0) = Mob*K->operator()(1,0);
	Kp(1,1) = Mob*K->operator()(1,1);
	Kp(1,2) = Mob*K->operator()(1,2);
    Kp(2,0) = Mob*K->operator()(2,0);
    Kp(2,1) = Mob*K->operator()(2,1);
    Kp(2,2) = Mob*K->operator()(2,2);

	// Loop on facets to build Rp, Np
	for(UInt localFacetId=0; localFacetId<numCellFacets; ++localFacetId)
    {
		const UInt globalFacetId = cellFacetsId[localFacetId];
        const Rigid_Mesh::Facet & fac = pMesh.getFacetsVector()[globalFacetId];

		const Point3D facetNormal     = fac.getUnsignedNormal();
        Point3D g                     = fac.getCentroid() - cellp.getCentroid();
		const Real facetMeasure       = fac.area();
		const Real alpha              = cellp.getAlpha(globalFacetId);
            
		Np(localFacetId,0) = facetNormal[0];
		Np(localFacetId,1) = facetNormal[1];
		Np(localFacetId,2) = facetNormal[2];
		Np *= Kp;

		g[0] = std::fabs(g[0]) < 1e-15 ? 0 : g[0];          // set relative eps
		g[1] = std::fabs(g[1]) < 1e-15 ? 0 : g[1];          // set relative eps
		g[2] = std::fabs(g[2]) < 1e-15 ? 0 : g[2];          // set relative eps

		Rp(localFacetId,0) = alpha * g[0] * facetMeasure;
		Rp(localFacetId,1) = alpha * g[1] * facetMeasure;
		Rp(localFacetId,2) = alpha * g[2] * facetMeasure;
	}

	// Building the Mp matrix
	Sp = Np.fullPivLu().image(Np);   // Because otherwise NtN could be singular in the unlucky case of 2 parallelel faces

	Mat     NtN = Sp.transpose() * Sp;
	Mat     NNtNiNt = Sp * NtN.inverse() * Sp.transpose();

	Mp0.resize(numCellFacets,numCellFacets);
	Mp0.setZero();
	Mp1.resize(numCellFacets,numCellFacets);
	Mp1.setZero();
	Mp.resize(numCellFacets,numCellFacets);
	Mp.setZero();

	Mp0 = ( 1. / cellMeasure ) *
		( Rp * Kp.inverse() * Rp.transpose() );
	Mp1 = Mp0.trace() * gamma / numCellFacets *
		( Eigen::MatrixXd::Identity(numCellFacets,numCellFacets) - NNtNiNt );

	Mp  = Mp0 + Mp1;
}

void local_InnerProduct::assemble_inv()      
{	
	Mat      Qp;         // I will take Np as the column base of Sp
	Mat33    Kp;         // permability tensor
	
	auto & K = pMesh.getPropertiesMap().getProperties(cellp.getZoneCode()).M_permeability; 
	Real Mob = pMesh.getPropertiesMap().getMobility();
	const std::vector<UInt> & cellFacetsId( cellp.getFacetsIds() );
	const UInt numCellFacets  = cellFacetsId.size();
    const Real cellMeasure    = cellp.getVolume();

    // Resize local matrices
    Np.resize(numCellFacets,Eigen::NoChange);
	Np.setZero();
    Rp.resize(numCellFacets,Eigen::NoChange);
	Rp.setZero();

	Kp(0,0) = Mob*K->operator()(0,0);
	Kp(0,1) = Mob*K->operator()(0,1);
	Kp(0,2) = Mob*K->operator()(0,2);
	Kp(1,0) = Mob*K->operator()(1,0);
	Kp(1,1) = Mob*K->operator()(1,1);
	Kp(1,2) = Mob*K->operator()(1,2);
    Kp(2,0) = Mob*K->operator()(2,0);
    Kp(2,1) = Mob*K->operator()(2,1);
    Kp(2,2) = Mob*K->operator()(2,2);

	// Loop on facets to build Rp, Np
	for(UInt localFacetId=0; localFacetId<numCellFacets; ++localFacetId)
    {
		const UInt globalFacetId = cellFacetsId[localFacetId];
        const Rigid_Mesh::Facet & fac = pMesh.getFacetsVector()[globalFacetId];

		const Point3D facetNormal     = fac.getUnsignedNormal();
        Point3D g                     = fac.getCentroid() - cellp.getCentroid();
		const Real facetMeasure       = fac.area();
		const Real alpha              = cellp.getAlpha(globalFacetId);
            
		Np(localFacetId,0) = facetNormal[0];
		Np(localFacetId,1) = facetNormal[1];
		Np(localFacetId,2) = facetNormal[2];
		Np *= Kp;

		g[0] = std::fabs(g[0]) < 1e-15 ? 0 : g[0];          // set relative eps
		g[1] = std::fabs(g[1]) < 1e-15 ? 0 : g[1];          // set relative eps
		g[2] = std::fabs(g[2]) < 1e-15 ? 0 : g[2];          // set relative eps

		Rp(localFacetId,0) = alpha * g[0] * facetMeasure;
		Rp(localFacetId,1) = alpha * g[1] * facetMeasure;
		Rp(localFacetId,2) = alpha * g[2] * facetMeasure;
	}

	// Building the Mp matrix
	Qp = Rp.fullPivLu().image(Rp);   // Because otherwise RtR could be singular in the unlucky case of 2 parallelel faces

    Mat RtR = Qp.transpose() * Qp;
    Mat RRtRiRt = Qp * RtR.inverse() * Qp.transpose();

	Mp0.resize(numCellFacets,numCellFacets);
	Mp0.setZero();
	Mp1.resize(numCellFacets,numCellFacets);
	Mp1.setZero();
	Mp.resize(numCellFacets,numCellFacets);
	Mp.setZero();

	Mp0 = ( 1. / cellMeasure ) *
		(Np * Kp.inverse() * Np.transpose());
	Mp1 = Mp0.trace() * gamma / numCellFacets *
		( Eigen::MatrixXd::Identity(numCellFacets,numCellFacets) - RRtRiRt );

    Mp  = Mp0 + Mp1;
}	


void local_Div::assemble()
{	
	const std::vector<UInt> & cellFacetsId( cellp.getFacetsIds() );
	const UInt numCellFacets  = cellFacetsId.size();

	// Loop on facets to build Rp, Np
	for(UInt localFacetId=0; localFacetId<numCellFacets; ++localFacetId)
    {
		const UInt globalFacetId = cellFacetsId[localFacetId];
        const Rigid_Mesh::Facet & fac = pMesh.getFacetsVector()[globalFacetId];
		const Real alpha              = cellp.getAlpha(globalFacetId);
		const Real facetMeasure       = fac.area();

		Bp.push_back( - alpha * facetMeasure );    // The minus because I have changed the sign to the conservation equation
	}
}


void local_builder::build()
{	
	Mat     Sp;          // I will take Np as the column base of Sp
	Mat33   Kp;          // permability tensor
	
	auto & K = IP.pMesh.getPropertiesMap().getProperties( cel.getZoneCode() ).M_permeability;
	const std::vector<UInt> & cellFacetsId( cel.getFacetsIds() );
	const UInt numCellFacets  = cellFacetsId.size();
    const Real cellMeasure    = cel.getVolume();

    // Resize local matrices
    IP.Np.resize(numCellFacets,Eigen::NoChange);
	IP.Np.setZero();
    IP.Rp.resize(numCellFacets,Eigen::NoChange);
	IP.Rp.setZero();

	Kp(0,0) = K->operator()(0,0);
	Kp(0,1) = K->operator()(0,1);
	Kp(0,2) = K->operator()(0,2);
	Kp(1,0) = K->operator()(1,0);
	Kp(1,1) = K->operator()(1,1);
	Kp(1,2) = K->operator()(1,2);
    Kp(2,0) = K->operator()(2,0);
    Kp(2,1) = K->operator()(2,1);
    Kp(2,2) = K->operator()(2,2);

	// Loop on facets to build Rp, Np
	for(UInt localFacetId=0; localFacetId<numCellFacets; ++localFacetId)
    {
		const UInt globalFacetId = cellFacetsId[localFacetId];
        const Rigid_Mesh::Facet & fac = IP.pMesh.getFacetsVector()[globalFacetId];

		const Point3D facetNormal     = fac.getUnsignedNormal();
        Point3D g                     = fac.getCentroid() - cel.getCentroid();
		const Real facetMeasure       = fac.area();
		const Real alpha              = cel.getAlpha(globalFacetId);
            
		IP.Np(localFacetId,0) = facetNormal[0];
		IP.Np(localFacetId,1) = facetNormal[1];
		IP.Np(localFacetId,2) = facetNormal[2];
		IP.Np *= Kp;

		g[0] = std::fabs(g[0]) < 1e-15 ? 0 : g[0];          // set relative eps
		g[1] = std::fabs(g[1]) < 1e-15 ? 0 : g[1];          // set relative eps
		g[2] = std::fabs(g[2]) < 1e-15 ? 0 : g[2];          // set relative eps

		IP.Rp(localFacetId,0) = alpha * g[0] * facetMeasure;
		IP.Rp(localFacetId,1) = alpha * g[1] * facetMeasure;
		IP.Rp(localFacetId,2) = alpha * g[2] * facetMeasure;
		
		Div.Bp.push_back( - alpha * facetMeasure );    // The minus because I have changed the sign to the conservation equation
	}

	// Building the Mp matrix
	Sp = IP.Np.fullPivLu().image(IP.Np);    // Because otherwise NtN could be singular in the unlucky case of 2 parallelel faces
//	Sp = IP.Np;

	Mat      NtN = Sp.transpose() * Sp;
	Mat      NNtNiNt = Sp * NtN.inverse() * Sp.transpose();

	IP.Mp0.resize(numCellFacets,numCellFacets);
	IP.Mp0.setZero();
	IP.Mp1.resize(numCellFacets,numCellFacets);
	IP.Mp1.setZero();
	IP.Mp.resize(numCellFacets,numCellFacets);
	IP.Mp.setZero();

	IP.Mp0 = ( 1. / cellMeasure ) *
		( IP.Rp * Kp.inverse() * IP.Rp.transpose() );
	IP.Mp1 = IP.Mp0.trace() * IP.gamma / numCellFacets *
		( Eigen::MatrixXd::Identity(numCellFacets,numCellFacets) - NNtNiNt );
//	IP.Mp1 = IP.Mp0.trace() * IP.gamma *
//		( Eigen::MatrixXd::Identity(numCellFacets,numCellFacets) - NNtNiNt );

	IP.Mp  = IP.Mp0 + IP.Mp1;

}	

}
