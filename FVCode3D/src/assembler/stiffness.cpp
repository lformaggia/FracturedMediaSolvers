/*!
 *	@file stiffness.cpp
 *	@brief This class build a Stiffness-matrix of the Darcy problem (definitions).
 */

#include "assembler/stiffness.hpp"
#include "property/Properties.hpp"
#include "Eigen/Core"
#include "Eigen/SparseCore"
#include "geometry/Point3D.hpp"
#include <iostream>
namespace Darcy
{

StiffMatrix::Generic_Point StiffMatrix::border_center(Fracture_Juncture fj) const
{
	return (this->M_mesh.getNodesVector()[fj.first] +
				(this->M_mesh.getNodesVector()[fj.second] -
				 this->M_mesh.getNodesVector()[fj.first]
				)/2.
		   );
}

Real StiffMatrix::Findfracturesalpha (const Fracture_Juncture fj, const UInt n_Id) const
{
	Generic_Point borderCenter = border_center(fj);
	Generic_Point cellCenter = this->M_mesh.getFractureFacetsIdsVector()[n_Id].getCentroid();

	Real A = M_properties.getProperties(this->M_mesh.getFractureFacetsIdsVector()[n_Id].getZoneCode()).M_aperture * (this->M_mesh.getNodesVector()[fj.second]-this->M_mesh.getNodesVector()[fj.first]).norm();
	Real k = M_properties.getProperties(this->M_mesh.getFractureFacetsIdsVector()[n_Id].getZoneCode()).M_permeability;
	Generic_Vector f;
	Real alpha;
	Real D;
	f = borderCenter - cellCenter;
	D = sqrt(dotProduct(f, f));
	f /= D;
	Generic_Vector normal = crossProduct(f, this->M_mesh.getNodesVector()[fj.second]-this->M_mesh.getNodesVector()[fj.first]); // k = f x l
	normal = crossProduct(this->M_mesh.getNodesVector()[fj.second]-this->M_mesh.getNodesVector()[fj.first], normal); // n = l x k
	normal.normalize();
	Real scalprod = fabs(dotProduct(normal, f));
	alpha = A*k*scalprod/D;
	return alpha;
}

Real StiffMatrix::Findalpha (const UInt & cellId, Facet_ID * const facet) const
{
	Generic_Point facetCenter = facet->getCentroid();
	Generic_Point cellCenter = this->M_mesh.getCellsVector()[cellId].getCentroid();
 	Generic_Vector normal = facet->getUNormal();
	Real alpha;
	Real D;
	Generic_Vector f;
	Real scalprod;
	f = cellCenter - facetCenter;
	D = sqrt(dotProduct(f, f));
	f = f/D;

	scalprod = fabs(dotProduct(f, normal));
	alpha = facet->getSize() * M_properties.getProperties(this->M_mesh.getCellsVector()[cellId].getZoneCode()).M_permeability * scalprod / D;
	return alpha;
}

Real StiffMatrix::FindDirichletalpha (const UInt & cellId, Facet_ID * const facet) const
{
	Generic_Point facetCenter = facet->getCentroid();
	Generic_Point cellCenter = this->M_mesh.getCellsVector()[cellId].getCentroid();
	Real alpha;
	Real D;
	Generic_Vector f;
	f = cellCenter - facetCenter;
	D = sqrt(dotProduct(f, f));

	alpha = facet->getSize() * M_properties.getProperties(this->M_mesh.getCellsVector()[cellId].getZoneCode()).M_permeability / D;
	return alpha;
}

/*
void StiffMatrix::assemble()
{
	std::vector<Triplet> Matrix_elements;
	Real alpha1, alpha2, alphaF, alphaf;
	UInt neighbor1id, neighbor2id;
	std::vector<Real> alphas;
	Generic_Vector f1, f2;
	Real Df;
	Real T12, Q12, T1f, Q1f, T2f, Q2f, QFf, Q1o;

	for (auto facet_it : this->M_mesh.getInternalFacetsIdsVector())
	{
		neighbor1id = facet_it.getSeparated()[0];
		neighbor2id = facet_it.getSeparated()[1];

		alpha1 = Findalpha(neighbor1id, &facet_it);
		alpha2 = Findalpha(neighbor2id, &facet_it);

		T12 = alpha1*alpha2/(alpha1 + alpha2);
		Q12 = T12 * M_properties.getMobility();

		Matrix_elements.emplace_back(Triplet (neighbor1id, neighbor1id, Q12));
		Matrix_elements.emplace_back(Triplet (neighbor1id, neighbor2id, -Q12));
		Matrix_elements.emplace_back(Triplet (neighbor2id, neighbor1id, -Q12));
		Matrix_elements.emplace_back(Triplet (neighbor2id, neighbor2id, Q12));
	}

	for (auto facet_it : this->M_mesh.getFractureFacetsIdsVector())
	{
		Real F_permeability = M_properties.getProperties(facet_it.getZoneCode()).M_permeability;
		Real F_aperture = M_properties.getProperties(facet_it.getZoneCode()).M_aperture;
		Df = F_aperture/2.;
		alphaf = facet_it.getSize() * F_permeability / Df;

		neighbor1id = facet_it.getSeparated()[0];
		neighbor2id = facet_it.getSeparated()[1];

		alpha1 = Findalpha(neighbor1id, &facet_it);
		alpha2 = Findalpha(neighbor2id, &facet_it);

		T1f = alpha1*alphaf/(alpha1 + alphaf);
		Q1f = T1f * M_properties.getMobility();

		T2f = alphaf*alpha2/(alphaf + alpha2);
		Q2f = T2f * M_properties.getMobility();

		Matrix_elements.emplace_back(Triplet (neighbor1id, neighbor1id, Q1f));
		Matrix_elements.emplace_back(Triplet (neighbor1id, facet_it.getIdasCell(), -Q1f));
		Matrix_elements.emplace_back(Triplet (facet_it.getIdasCell(), neighbor1id, -Q1f));
		Matrix_elements.emplace_back(Triplet (facet_it.getIdasCell(), facet_it.getIdasCell(), Q1f));

		Matrix_elements.emplace_back(Triplet (neighbor2id, neighbor2id, Q2f));
		Matrix_elements.emplace_back(Triplet (neighbor2id, facet_it.getIdasCell(), -Q2f));
		Matrix_elements.emplace_back(Triplet (facet_it.getIdasCell(), neighbor2id, -Q2f));
		Matrix_elements.emplace_back(Triplet (facet_it.getIdasCell(), facet_it.getIdasCell(), Q2f));

		for (auto juncture_it : facet_it.getFractureNeighbors())
		{
			alphaF = Findfracturesalpha (juncture_it.first, facet_it.getId());

			for (auto neighbors_it : juncture_it.second)
				alphas.emplace_back(Findfracturesalpha (juncture_it.first, neighbors_it));

			Real a_sum = alphaF;
			auto sum_maker = [&a_sum](Real item){a_sum += item;};
			std::for_each(alphas.begin(), alphas.end(), sum_maker);

			for (UInt counter = 0; counter < alphas.size(); ++counter)
			{
				QFf = alphaF*alphas[counter] * M_properties.getMobility() / a_sum;
				Matrix_elements.emplace_back(Triplet (facet_it.getIdasCell(), facet_it.getIdasCell(), QFf));
				Matrix_elements.emplace_back(Triplet (facet_it.getIdasCell(), this->M_mesh.getFractureFacetsIdsVector()[juncture_it.second[counter]].getIdasCell(), -QFf));
			}
			alphas.clear();
		}
	}

	for(UInt i=0; i<this->M_size; ++i)
		_b->operator()(i)=0.;

	for (auto facet_it : this->M_mesh.getBorderFacetsIdsVector())
	{
		neighbor1id = facet_it.getSeparated()[0];
		UInt borderId = facet_it.getBorderId();

		if(m_Bc.getBordersBCMap().at(borderId).getBCType() == Neumann )
		{
			// Nella cond di bordo c'è già il contributo della permeabilità e mobilità (ma non la densità!)
			Q1o = m_Bc.getBordersBCMap().at(borderId).getBC()(facet_it.getCentroid()) * facet_it.getSize();

			_b->operator()(neighbor1id) += Q1o;
		}
		else if(m_Bc.getBordersBCMap().at(borderId).getBCType() == Dirichlet)
		{
			alpha1 = Findalpha (neighbor1id, &facet_it);
			alpha2 = FindDirichletalpha (neighbor1id, &facet_it);

			T12 = 2 * alpha1*alpha2/(alpha1 + alpha2);
			Q12 = T12 * M_properties.getMobility();
			Q1o = Q12 * m_Bc.getBordersBCMap().at(borderId).getBC()(facet_it.getCentroid());

			Matrix_elements.emplace_back(Triplet (neighbor1id, neighbor1id, Q12));
			_b->operator()(neighbor1id) = _b->operator()(neighbor1id) + Q1o;

		}
	}

	this->M_Matrix->setFromTriplets(Matrix_elements.begin(), Matrix_elements.end());
}
*/
//void StiffMatrix::assembleMFD(Real tCoeff)
void StiffMatrix::assemble()
{
	Real tCoeff=6.;
	using Eigen::Dynamic;
	using Eigen::RowMajor;
	using Eigen::ColMajor;
	using Geometry::Rigid_Mesh;
	std::vector<Triplet> Matrix_elements;
	std::vector<UInt> ZMatrix_elements;
	std::vector<UInt> BMatrix_elements;
	Real alpha1, alpha2, alphaF, alphaf;
	UInt neighbor1id, neighbor2id;
	std::vector<Real> alphas;
	Generic_Vector f1, f2;
	Real Df;
	Real T12, Q12, T1f, Q1f, T2f, Q2f, QFf, Q1o;
	// Local Matrices for Mimetic
	Eigen::Matrix<double,Dynamic,2> Np; // facet normals * K
	Eigen::Matrix<double,Dynamic,2> Nop; // facets normals outward
	Eigen::Matrix<double,Dynamic,Dynamic> Z0p;// COmpunent of Z matrix
	Eigen::Matrix<double,Dynamic,Dynamic> Z1p;// Component of Z Matrix
	Eigen::Matrix<double,Dynamic,Dynamic> Zp; //! Z matrix for internal product= \f$ M^{-1}\f$
	Eigen::Matrix<double,Dynamic,Dynamic> Tp; // Local trasmissibility
	Eigen::Matrix<double,1,Dynamic> Bp; // B Matrix
	Eigen::Matrix<double,Dynamic,1> BTp;// Traspose of B
	Eigen::Matrix<double,Dynamic,2> Rp; // Matrix R
	Eigen::Matrix<double,2,2> Rrt;
	Eigen::Matrix<double,2,2> Rinv;

	// Extract info of mesh. auto&& resolves const &
	auto&& cellVectorRef   = this->M_mesh.getCellsVector();
	auto&& facetVectorRef  = this->M_mesh.getFacetsVector();
	UInt numFacets = facetVectorRef.size();
//	UInt numInternalFacets=this->M_mesh.getInternalFacetsIdsVector().size();
	UInt numCells = cellVectorRef.size();

	// SIZING GLOBAL MATRICES
	Eigen::SparseMatrix<double,RowMajor> Z(numFacets,numFacets); //! Z matrix for internal product= \f$ M^{-1}\f$
	Eigen::SparseMatrix<double,RowMajor> B(numCells,numFacets);
	Eigen::SparseMatrix<double,RowMajor> T(numCells,numCells);
	ZMatrix_elements.resize(numFacets,0);
	BMatrix_elements.resize(numCells,0);
	// First loop to size matrices and avoid memory realloc
	for (auto&& cell : cellVectorRef)
	{
		std::vector<UInt> const & cellFacetsId( cell.getFacetsIds() );
		UInt numCellFacets = cellFacetsId.size();
		auto cellId         = cell.getId();
#ifndef NDEBUG
		if (cellId>numCells) std::cerr<<" INCORRECT CELL ID"<<std::endl;
#endif

		for(UInt localFacetId=0; localFacetId<numCellFacets;++localFacetId)
		{
			UInt globalFacetId = cellFacetsId[localFacetId];
			Rigid_Mesh::Facet const & fac=facetVectorRef[globalFacetId];
#ifndef NDEBUG
			if (globalFacetId>numFacets) std::cerr<<" INCORRECT FACE ID"<<std::endl;
#endif
			if(fac.isBorderFacet()) continue;
			BMatrix_elements[cellId]+=1;
#ifdef DIAGONALZ
			ZMatrix_elements[globalFacetId]+=1;
#else
			ZMatrix_elements[globalFacetId]+=numCellFacets;
#endif
			}
	}
	Z.reserve(ZMatrix_elements);
	B.reserve(BMatrix_elements);
	ZMatrix_elements.clear();
	ZMatrix_elements.shrink_to_fit();
	BMatrix_elements.clear();
	BMatrix_elements.shrink_to_fit();
	// Loop on cells
	std::cout<<" Starting loops on cells"<<std::endl;
	for (auto&& cell : cellVectorRef)
	{
		std::vector<UInt> const & cellFacetsId( cell.getFacetsIds() );
		UInt numCellFacets = cellFacetsId.size();
		auto K = M_properties.getMobility();
		auto cellBaricenter = cell.getCentroid();
		auto cellMeasure    = cell.getVolume();
		auto cellId         = cell.getId();
		if(cellId % 500 == 0)std::cout<<"Done "<< cellId<<" Cells"<<std::endl;
		// Resize local matrices
		Np.resize(numCellFacets,Eigen::NoChange);
		Np.setZero();
		Rp.resize(numCellFacets,Eigen::NoChange);
		Rp.setZero();
		Zp.resize(numCellFacets,numCellFacets);
		Zp.setZero();
		Z0p.resize(numCellFacets,numCellFacets);
		Z0p.setZero();
		Z1p.resize(numCellFacets,numCellFacets);
		Z1p.setZero();
		Bp.resize(Eigen::NoChange,numCellFacets);
		Bp.setZero();
		for(UInt localFacetId=0; localFacetId<numCellFacets;++localFacetId)
		{
			UInt globalFacetId = cellFacetsId[localFacetId];
			Rigid_Mesh::Facet const & fac=facetVectorRef[globalFacetId];
			if(fac.isBorderFacet()) continue;
			Real alpha(0.);
			Geometry::Point3D facetBaricenter= fac.getCentroid();
			Geometry::Point3D facetNormal    = fac.getUnsignedNormal();
			auto facetMeasure   = fac.area();
			Geometry::Point3D g = facetBaricenter-cellBaricenter;
			Real  dotp= Geometry::dotProduct(g,facetNormal);
			alpha = (dotp >=0.? 1.0:-1.0);
			//! @todo I use the formulation of Nicola (to be reviewed)
			// BEWARE FOR THE CONDITION ON VELOCITY I NEED THE AVERAGE
			// VELOCITY NOT THE FLUX
			Np(localFacetId,0)=facetNormal[0];
			Np(localFacetId,1)=facetNormal[1];
			Rp(localFacetId,0)=alpha*g[0]*facetMeasure;
			Rp(localFacetId,1)=alpha*g[1]*facetMeasure;
			Bp(0,localFacetId)=alpha*facetMeasure;
		}
		Rrt = Rp.transpose()*Rp;
		Real det  = 1.0/(Rrt(0,0)*Rrt(1,1)-2*Rrt(0,1));
		Rinv(0,0) =  det*Rrt(1,1);
		Rinv(0,1) = -det*Rrt(0,1);
		Rinv(1,0) =  Rinv(0,1);
		Rinv(1,1) =  det*Rrt(0,0);

		Z0p = (K/cellMeasure)*(Np*Np.transpose());
		Z1p = (tCoeff*K/cellMeasure)*(
				Eigen::MatrixXd::Identity(numCellFacets,numCellFacets)-
				Rp*Rinv*Rp.transpose());
		Zp  = Z0p + Z1p;
		for(UInt iloc=0; iloc<numCellFacets;++iloc)
		{
			UInt i=cellFacetsId[iloc];
			Real Bcoeff=Bp(0,iloc);
			if(Bcoeff != 0.) B.insert(cellId,i)=Bcoeff;
			Real Zcoeff = Zp(iloc,iloc);
			if(Zcoeff != 0. ) Z.coeffRef(i,i)+=Zp(iloc,iloc);
			for(UInt jloc=iloc+1; jloc<numCellFacets;++jloc)
			{
				UInt j=cellFacetsId[jloc];
				Zcoeff = Zp(iloc,jloc);
				if (Zcoeff != 0.)
				{
#ifdef DIAGONALZ
					Z.coeffRef(i,i)+=Zcoeff;
					Z.coeffRef(j,j)+=Zcoeff;
#else
					Z.coeffRef(i,j)+=Zcoeff;
					Z.coeffRef(j,i)+=Zcoeff;
#endif
				}
			}
		}
	}
// Now I Need to build T
	std::cout<<"Building Trasmissibility matrix (global)"<<std::endl;
	T = B * Z * B.transpose();
	T.makeCompressed();
// Clean up to seve memory (I hope)
	std::cout<<" Matrix T has "<<T.rows()<<" rows and "<<T.cols()<<" columns and "<<T.nonZeros()<<" Non zeros"<<std::endl;
	std::cout<<" Num Cells"<<numCells<<std::endl;
	std::cout.flush();
	B.resize(0,0); // Save memory but it does not work...
	Z.resize(0,0);

	std::cout<<"Building the triplets...."; std::cout.flush();
	UInt counter(0);
	Matrix_elements.reserve(T.nonZeros()+8*this->M_mesh.getFractureFacetsIdsVector().size());
	for (int k=0; k<T.outerSize(); ++k)
	for (Eigen::SparseMatrix<double>::InnerIterator it(T,k); it; ++it)
	{
		if (++counter % 1000 ==0) std::cerr<<"Doing: "<<counter<<std::endl;std::cout.flush();
		Matrix_elements.emplace_back(
				Triplet (it.row(), it.col(), it.value()));
	}
	T.resize(0,0);
	std::cout<<" Done ALL"<<std::endl;
	std::cout.flush();
	for (auto facet_it : this->M_mesh.getFractureFacetsIdsVector())
	{
		Real F_permeability = M_properties.getProperties(facet_it.getZoneCode()).M_permeability;
		Real F_aperture = M_properties.getProperties(facet_it.getZoneCode()).M_aperture;
		Df = F_aperture/2.;
		alphaf = facet_it.getSize() * F_permeability / Df;

		neighbor1id = facet_it.getSeparated()[0];
		neighbor2id = facet_it.getSeparated()[1];

		alpha1 = Findalpha(neighbor1id, &facet_it);
		alpha2 = Findalpha(neighbor2id, &facet_it);

		T1f = alpha1*alphaf/(alpha1 + alphaf);
		Q1f = T1f * M_properties.getMobility();

		T2f = alphaf*alpha2/(alphaf + alpha2);
		Q2f = T2f * M_properties.getMobility();

		Matrix_elements.emplace_back(Triplet (neighbor1id, neighbor1id, Q1f));
		Matrix_elements.emplace_back(Triplet (neighbor1id, facet_it.getIdasCell(), -Q1f));
		Matrix_elements.emplace_back(Triplet (facet_it.getIdasCell(), neighbor1id, -Q1f));
		Matrix_elements.emplace_back(Triplet (facet_it.getIdasCell(), facet_it.getIdasCell(), Q1f));

		Matrix_elements.emplace_back(Triplet (neighbor2id, neighbor2id, Q2f));
		Matrix_elements.emplace_back(Triplet (neighbor2id, facet_it.getIdasCell(), -Q2f));
		Matrix_elements.emplace_back(Triplet (facet_it.getIdasCell(), neighbor2id, -Q2f));
		Matrix_elements.emplace_back(Triplet (facet_it.getIdasCell(), facet_it.getIdasCell(), Q2f));

		for (auto juncture_it : facet_it.getFractureNeighbors())
		{
			alphaF = Findfracturesalpha (juncture_it.first, facet_it.getId());

			for (auto neighbors_it : juncture_it.second)
				alphas.emplace_back(Findfracturesalpha (juncture_it.first, neighbors_it));

			Real a_sum = alphaF;
			auto sum_maker = [&a_sum](Real item){a_sum += item;};
			std::for_each(alphas.begin(), alphas.end(), sum_maker);

			for (UInt counter = 0; counter < alphas.size(); ++counter)
			{
				QFf = alphaF*alphas[counter] * M_properties.getMobility() / a_sum;
				Matrix_elements.emplace_back(Triplet (facet_it.getIdasCell(), facet_it.getIdasCell(), QFf));
				Matrix_elements.emplace_back(Triplet (facet_it.getIdasCell(), this->M_mesh.getFractureFacetsIdsVector()[juncture_it.second[counter]].getIdasCell(), -QFf));
			}
			alphas.clear();
		}
	}

	for(UInt i=0; i<this->M_size; ++i)
		_b->operator()(i)=0.;

	for (auto facet_it : this->M_mesh.getBorderFacetsIdsVector())
	{
		neighbor1id = facet_it.getSeparated()[0];
		UInt borderId = facet_it.getBorderId();

		if(m_Bc.getBordersBCMap().at(borderId).getBCType() == Neumann )
		{
			// Nella cond di bordo c'è già il contributo della permeabilità e mobilità (ma non la densità!)
			Q1o = m_Bc.getBordersBCMap().at(borderId).getBC()(facet_it.getCentroid()) * facet_it.getSize();

			_b->operator()(neighbor1id) += Q1o;
		}
		else if(m_Bc.getBordersBCMap().at(borderId).getBCType() == Dirichlet)
		{
			alpha1 = Findalpha (neighbor1id, &facet_it);
			alpha2 = FindDirichletalpha (neighbor1id, &facet_it);

			T12 = 2 * alpha1*alpha2/(alpha1 + alpha2);
			Q12 = T12 * M_properties.getMobility();
			Q1o = Q12 * m_Bc.getBordersBCMap().at(borderId).getBC()(facet_it.getCentroid());

			Matrix_elements.emplace_back(Triplet (neighbor1id, neighbor1id, Q12));
			_b->operator()(neighbor1id) = _b->operator()(neighbor1id) + Q1o;

		}
	}

	this->M_Matrix->setFromTriplets(Matrix_elements.begin(), Matrix_elements.end());
	std::cout<<" Assembling ended"<<std::endl;
}

} // namespace Darcy
