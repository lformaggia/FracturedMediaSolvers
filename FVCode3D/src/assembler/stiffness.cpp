/*!
 *	@file stiffness.cpp
 *	@brief This class build a Stiffness-matrix of the Darcy problem (definitions).
 */

#include "assembler/stiffness.hpp"
#include "property/Properties.hpp"

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

			T12 = alpha1*alpha2/(alpha1 + alpha2);
			Q12 = T12 * M_properties.getMobility();
			Q1o = Q12 * m_Bc.getBordersBCMap().at(borderId).getBC()(facet_it.getCentroid());

			Matrix_elements.emplace_back(Triplet (neighbor1id, neighbor1id, Q12));
			_b->operator()(neighbor1id) = _b->operator()(neighbor1id) + Q1o;

		}
	}

	this->M_Matrix->setFromTriplets(Matrix_elements.begin(), Matrix_elements.end());
}

} // namespace Darcy
