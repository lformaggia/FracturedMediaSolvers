/*!
 *  @file stiffness.cpp
 *  @brief This class build a Stiffness-matrix of the Darcy problem (definitions).
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

    if( !this->M_mesh.getInternalFacetsIdsVector().empty() )
    {
        assemblePorousMatrix( Matrix_elements );
    } // if

    if( !this->M_mesh.getFractureFacetsIdsVector().empty() )
    {
        assembleFracture( Matrix_elements );
    } // if

    if( !this->M_mesh.getInternalFacetsIdsVector().empty() )
    {
        assemblePorousMatrixBC( Matrix_elements );
    } // if

    this->M_Matrix->setFromTriplets( std::begin( Matrix_elements ), std::end( Matrix_elements ) );
} // assemble

void StiffMatrix::assemblePorousMatrix( std::vector<Triplet>& Matrix_elements ) const
{
    for (auto facet_it : this->M_mesh.getInternalFacetsIdsVector())
    {
        const UInt neighbor1id = facet_it.getSeparated()[0];
        const UInt neighbor2id = facet_it.getSeparated()[1];

        const Real alpha1 = Findalpha(neighbor1id, &facet_it);
        const Real alpha2 = Findalpha(neighbor2id, &facet_it);

        const Real T12 = alpha1*alpha2/(alpha1 + alpha2);
        const Real Q12 = T12 * M_properties.getMobility();

        Matrix_elements.emplace_back(Triplet (neighbor1id, neighbor1id, Q12));
        Matrix_elements.emplace_back(Triplet (neighbor1id, neighbor2id, -Q12));
        Matrix_elements.emplace_back(Triplet (neighbor2id, neighbor1id, -Q12));
        Matrix_elements.emplace_back(Triplet (neighbor2id, neighbor2id, Q12));
    } // for

    for (auto facet_it : this->M_mesh.getFractureFacetsIdsVector())
    {
        const Real F_permeability = M_properties.getProperties(facet_it.getZoneCode()).M_permeability;
        const Real F_aperture = M_properties.getProperties(facet_it.getZoneCode()).M_aperture;
        const Real Df = F_aperture/2.;
        const Real alphaf = facet_it.getSize() * F_permeability / Df;

        const UInt neighbor1id = facet_it.getSeparated()[0];
        const UInt neighbor2id = facet_it.getSeparated()[1];

        const Real alpha1 = Findalpha(neighbor1id, &facet_it);
        const Real alpha2 = Findalpha(neighbor2id, &facet_it);

        const Real T1f = alpha1*alphaf/(alpha1 + alphaf);
        const Real Q1f = T1f * M_properties.getMobility();

        const Real T2f = alphaf*alpha2/(alphaf + alpha2);
        const Real Q2f = T2f * M_properties.getMobility();

        Matrix_elements.emplace_back(Triplet (neighbor1id, neighbor1id, Q1f));
        Matrix_elements.emplace_back(Triplet (neighbor1id, facet_it.getIdasCell(), -Q1f));
        Matrix_elements.emplace_back(Triplet (facet_it.getIdasCell(), neighbor1id, -Q1f));
        Matrix_elements.emplace_back(Triplet (facet_it.getIdasCell(), facet_it.getIdasCell(), Q1f));

        Matrix_elements.emplace_back(Triplet (neighbor2id, neighbor2id, Q2f));
        Matrix_elements.emplace_back(Triplet (neighbor2id, facet_it.getIdasCell(), -Q2f));
        Matrix_elements.emplace_back(Triplet (facet_it.getIdasCell(), neighbor2id, -Q2f));
        Matrix_elements.emplace_back(Triplet (facet_it.getIdasCell(), facet_it.getIdasCell(), Q2f));
    } // for
} // StiffMatrix::assemblePorousMatrix

void StiffMatrix::assemblePorousMatrixBC( std::vector<Triplet>& Matrix_elements ) const
{
    for(UInt i=0; i<this->M_size; ++i)
    {
        _b->operator()(i)=0.;
    } // for

    for (auto facet_it : this->M_mesh.getBorderFacetsIdsVector())
    {
        const UInt neighbor1id = facet_it.getSeparated()[0];
        const UInt borderId = facet_it.getBorderId();

        if(m_Bc.getBordersBCMap().at(borderId).getBCType() == Neumann )
        {
            // Nella cond di bordo c'è già il contributo della permeabilità e mobilità (ma non la densità!)
            const Real Q1o = m_Bc.getBordersBCMap().at(borderId).getBC()(facet_it.getCentroid()) * facet_it.getSize();

            _b->operator()(neighbor1id) += Q1o;
        } // if
        else if(m_Bc.getBordersBCMap().at(borderId).getBCType() == Dirichlet)
        {
            const Real alpha1 = Findalpha (neighbor1id, &facet_it);
            const Real alpha2 = FindDirichletalpha (neighbor1id, &facet_it);

            const Real T12 = alpha1*alpha2/(alpha1 + alpha2);
            const Real Q12 = T12 * M_properties.getMobility();
            const Real Q1o = Q12 * m_Bc.getBordersBCMap().at(borderId).getBC()(facet_it.getCentroid());

            Matrix_elements.emplace_back(Triplet (neighbor1id, neighbor1id, Q12));
            _b->operator()(neighbor1id) = _b->operator()(neighbor1id) + Q1o;
        } // else if
    } // for
} // StiffMatrix::assemblePorousMatrixBC

void StiffMatrix::assembleFracture( std::vector<Triplet>& Matrix_elements ) const
{
    for (auto facet_it : this->M_mesh.getFractureFacetsIdsVector())
    {
        for (auto juncture_it : facet_it.getFractureNeighbors())
        {
            std::vector<Real> alphas;
            const Real alphaF = Findfracturesalpha (juncture_it.first, facet_it.getId());

            for (auto neighbors_it : juncture_it.second)
            {
                alphas.emplace_back(Findfracturesalpha (juncture_it.first, neighbors_it));
            } // for

            Real a_sum = alphaF;
            auto sum_maker = [&a_sum](Real item){a_sum += item;};
            std::for_each(alphas.begin(), alphas.end(), sum_maker);

            for (UInt counter = 0; counter < alphas.size(); ++counter)
            {
                const Real QFf = alphaF*alphas[counter] * M_properties.getMobility() / a_sum;
                Matrix_elements.emplace_back(Triplet (facet_it.getIdasCell(), facet_it.getIdasCell(), QFf));
                Matrix_elements.emplace_back(Triplet (facet_it.getIdasCell(),
                    this->M_mesh.getFractureFacetsIdsVector()[juncture_it.second[counter]].getIdasCell(), -QFf));
            } // for
        } // for
    } // for
} // StiffMatrix::assembleFracture

} // namespace Darcy
