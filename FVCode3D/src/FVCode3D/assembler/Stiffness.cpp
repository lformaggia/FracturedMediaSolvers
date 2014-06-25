/*!
 *  @file stiffness.cpp
 *  @brief This class build a Stiffness-matrix of the Darcy problem (definitions).
 */

#include <FVCode3D/assembler/Stiffness.hpp>
#include <FVCode3D/property/Properties.hpp>

namespace FVCode3D
{

Point3D StiffMatrix::getBorderCenter(Fracture_Juncture fj) const
{
    return (this->M_mesh.getNodesVector()[fj.first] +
                (this->M_mesh.getNodesVector()[fj.second] -
                 this->M_mesh.getNodesVector()[fj.first]
                )/2.
           );
}

Real StiffMatrix::findFracturesAlpha (const Fracture_Juncture fj, const UInt n_Id) const
{
	Point3D borderCenter = getBorderCenter(fj);
	Point3D cellCenter = this->M_mesh.getFractureFacetsIdsVector()[n_Id].getCentroid();

    Real A = M_properties.getProperties(this->M_mesh.getFractureFacetsIdsVector()[n_Id].getZoneCode()).M_aperture * (this->M_mesh.getNodesVector()[fj.second]-this->M_mesh.getNodesVector()[fj.first]).norm();
    Real k = M_properties.getProperties(this->M_mesh.getFractureFacetsIdsVector()[n_Id].getZoneCode()).M_permeability;
    Point3D f;
    Real alpha;
    Real D;
    f = borderCenter - cellCenter;
    D = sqrt(dotProduct(f, f));
    f /= D;
    Point3D normal = crossProduct(f, this->M_mesh.getNodesVector()[fj.second]-this->M_mesh.getNodesVector()[fj.first]); // k = f x l
    normal = crossProduct(this->M_mesh.getNodesVector()[fj.second]-this->M_mesh.getNodesVector()[fj.first], normal); // n = l x k
    normal.normalize();
    Real scalprod = fabs(dotProduct(normal, f));
    alpha = A * k * scalprod / D;
    return alpha;
}

Real StiffMatrix::findAlpha (const UInt & cellId, const Facet_ID * facet) const
{
    Point3D facetCenter = facet->getCentroid();
    Point3D cellCenter = this->M_mesh.getCellsVector()[cellId].getCentroid();
    Point3D normal = facet->getUnsignedNormal();
    Real alpha;
    Real D;
    Point3D f;
    Real scalprod;
    f = cellCenter - facetCenter;
    D = sqrt(dotProduct(f, f));
    f = f/D;

    scalprod = fabs(dotProduct(f, normal));
    alpha = facet->getSize() * M_properties.getProperties(this->M_mesh.getCellsVector()[cellId].getZoneCode()).M_permeability * scalprod / D;
    return alpha;
}

Real StiffMatrix::findAlpha (const UInt & facetId, const Edge_ID * edge) const
{
    Point3D edgeCenter = edge->getCentroid();
    Point3D facetCenter = this->M_mesh.getFacetsVector()[facetId].getCentroid();

    Real A = edge->getSize() * M_properties.getProperties(this->M_mesh.getFacetsVector()[facetId].getZoneCode()).M_aperture;
    Real k = M_properties.getProperties(this->M_mesh.getFacetsVector()[facetId].getZoneCode()).M_permeability;
    Point3D f;
    Real alpha;
    Real D;
    f = edgeCenter - facetCenter;
    D = sqrt(dotProduct(f, f));
    f = f/D;

    Point3D normal = crossProduct(	f,
    						this->M_mesh.getNodesVector()[edge->getEdge().getVerticesIds()[0]] -
    						this->M_mesh.getNodesVector()[edge->getEdge().getVerticesIds()[1]] ); // k = f x l
    normal = crossProduct(	this->M_mesh.getNodesVector()[edge->getEdge().getVerticesIds()[0]] -
    						this->M_mesh.getNodesVector()[edge->getEdge().getVerticesIds()[1]],
    						normal); // n = l x k
    normal.normalize();
    Real scalprod = fabs(dotProduct(f, normal));
    alpha = A * k * scalprod / D;
	return alpha;
}

Real StiffMatrix::findDirichletAlpha (const UInt & cellId, const Facet_ID * facet) const
{
    Point3D facetCenter = facet->getCentroid();
    Point3D cellCenter = this->M_mesh.getCellsVector()[cellId].getCentroid();
    Real alpha;
    Real D;
    Point3D f = cellCenter - facetCenter;
    D = sqrt(dotProduct(f, f));

    alpha = facet->getSize() * M_properties.getProperties(this->M_mesh.getCellsVector()[cellId].getZoneCode()).M_permeability / D;
    return alpha;
}

Real StiffMatrix::findDirichletAlpha (const UInt & facetId, const Edge_ID * edge) const
{
    Point3D borderCenter = edge->getCentroid();
    Point3D facetCenter = this->M_mesh.getFacetsVector()[facetId].getCentroid();

    Real A = edge->getSize() * M_properties.getProperties(this->M_mesh.getFacetsVector()[facetId].getZoneCode()).M_aperture;
    Real k = M_properties.getProperties(this->M_mesh.getFacetsVector()[facetId].getZoneCode()).M_permeability;
    Real alpha;
    Real D;
    Point3D f = borderCenter - facetCenter;
    D = sqrt(dotProduct(f, f));

    alpha = A * k / D;
    return alpha;
}

void StiffMatrix::assemble()
{
    // TODO reserve the vector!

    if( !this->M_mesh.getInternalFacetsIdsVector().empty() )
    {
        assemblePorousMatrix();
    } // if

    if( !this->M_mesh.getFractureFacetsIdsVector().empty() )
    {
        assembleFracture();
    } // if

    if( !this->M_mesh.getInternalFacetsIdsVector().empty() )
    {
        assemblePorousMatrixBC();
    } // if

    if( !this->M_mesh.getBorderTipEdgesIdsVector().empty() )
    {
        assembleFractureBC();
    } // if
} // assemble

void StiffMatrix::assemblePorousMatrix()
{
    for (auto& facet_it : this->M_mesh.getInternalFacetsIdsVector())
    {
        const UInt neighbor1id = facet_it.getSeparatedCellsIds()[0];
        const UInt neighbor2id = facet_it.getSeparatedCellsIds()[1];

        const Real alpha1 = findAlpha(neighbor1id, &facet_it);
        const Real alpha2 = findAlpha(neighbor2id, &facet_it);

        const Real T12 = alpha1*alpha2/(alpha1 + alpha2);
        const Real Q12 = T12 * M_properties.getMobility();

        M_matrixElements.emplace_back(M_offsetRow + neighbor1id, M_offsetCol + neighbor1id, Q12); // Triplet
        M_matrixElements.emplace_back(M_offsetRow + neighbor1id, M_offsetCol + neighbor2id, -Q12); // Triplet
        M_matrixElements.emplace_back(M_offsetRow + neighbor2id, M_offsetCol + neighbor1id, -Q12); // Triplet
        M_matrixElements.emplace_back(M_offsetRow + neighbor2id, M_offsetCol + neighbor2id, Q12); // Triplet
    } // for

    for (auto& facet_it : this->M_mesh.getFractureFacetsIdsVector())
    {
        const Real F_permeability = M_properties.getProperties(facet_it.getZoneCode()).M_permeability;
        const Real F_aperture = M_properties.getProperties(facet_it.getZoneCode()).M_aperture;
        const Real Df = F_aperture/2.;
        const Real alphaf = facet_it.getSize() * F_permeability / Df;

        const UInt neighbor1id = facet_it.getSeparatedCellsIds()[0];
        const UInt neighbor2id = facet_it.getSeparatedCellsIds()[1];

        const Real alpha1 = findAlpha(neighbor1id, &facet_it);
        const Real alpha2 = findAlpha(neighbor2id, &facet_it);

        const Real T1f = alpha1*alphaf/(alpha1 + alphaf);
        const Real Q1f = T1f * M_properties.getMobility();

        const Real T2f = alphaf*alpha2/(alphaf + alpha2);
        const Real Q2f = T2f * M_properties.getMobility();

        M_matrixElements.emplace_back(M_offsetRow + neighbor1id, M_offsetCol + neighbor1id, Q1f); // Triplet
        M_matrixElements.emplace_back(M_offsetRow + neighbor1id, M_offsetCol + facet_it.getIdAsCell(), -Q1f); // Triplet
        M_matrixElements.emplace_back(M_offsetRow + facet_it.getIdAsCell(), M_offsetCol + neighbor1id, -Q1f); // Triplet
        M_matrixElements.emplace_back(M_offsetRow + facet_it.getIdAsCell(), M_offsetCol + facet_it.getIdAsCell(), Q1f); // Triplet

        M_matrixElements.emplace_back(M_offsetRow + neighbor2id, M_offsetCol + neighbor2id, Q2f); // Triplet
        M_matrixElements.emplace_back(M_offsetRow + neighbor2id, M_offsetCol + facet_it.getIdAsCell(), -Q2f); // Triplet
        M_matrixElements.emplace_back(M_offsetRow + facet_it.getIdAsCell(), M_offsetCol + neighbor2id, -Q2f); // Triplet
        M_matrixElements.emplace_back(M_offsetRow + facet_it.getIdAsCell(), M_offsetCol + facet_it.getIdAsCell(), Q2f); // Triplet
    } // for
} // StiffMatrix::assemblePorousMatrix

void StiffMatrix::assemblePorousMatrixBC()
{
    for (auto& facet_it : this->M_mesh.getBorderFacetsIdsVector())
    {
        const UInt neighbor1id = facet_it.getSeparatedCellsIds()[0];
        const UInt borderId = facet_it.getBorderId();

        if(M_bc.getBordersBCMap().at(borderId).getBCType() == Neumann )
        {
            // Nella cond di bordo c'è già il contributo della permeabilità e mobilità (ma non la densità!)
            const Real Q1o = M_bc.getBordersBCMap().at(borderId).getBC()(facet_it.getCentroid()) * facet_it.getSize();

            M_b->operator()(M_offsetRow + neighbor1id) += Q1o;
        } // if
        else if(M_bc.getBordersBCMap().at(borderId).getBCType() == Dirichlet)
        {
            const Real alpha1 = findAlpha (neighbor1id, &facet_it);
            const Real alpha2 = findDirichletAlpha (neighbor1id, &facet_it);

            const Real T12 = alpha1*alpha2/(alpha1 + alpha2);
            const Real Q12 = T12 * M_properties.getMobility();
            const Real Q1o = Q12 * M_bc.getBordersBCMap().at(borderId).getBC()(facet_it.getCentroid());

            M_matrixElements.emplace_back(M_offsetRow + neighbor1id, M_offsetCol + neighbor1id, Q12); // Triplet
            M_b->operator()(M_offsetRow + neighbor1id) += Q1o;
        } // else if
    } // for
} // StiffMatrix::assemblePorousMatrixBC

void StiffMatrix::assembleFracture()
{
    for(auto& facet_it : this->M_mesh.getFractureFacetsIdsVector())
    {
        for(auto& juncture_it : facet_it.getFractureNeighbors())
        {
            std::vector<Real> alphas;
            const Real alphaF = findFracturesAlpha(juncture_it.first, facet_it.getFractureId());

            for(auto neighbors_it : juncture_it.second)
            {
                alphas.emplace_back(findFracturesAlpha(juncture_it.first, neighbors_it));
            } // for

            Real a_sum = alphaF;
            auto sum_maker = [&a_sum](Real item){a_sum += item;};
            std::for_each(alphas.begin(), alphas.end(), sum_maker);

            for(UInt counter = 0; counter < alphas.size(); ++counter)
            {
                const Real QFf = alphaF*alphas[counter] * M_properties.getMobility() / a_sum;
                M_matrixElements.emplace_back(M_offsetRow + facet_it.getIdAsCell(), M_offsetCol + facet_it.getIdAsCell(), QFf); // Triplet
                M_matrixElements.emplace_back(M_offsetRow + facet_it.getIdAsCell(),
                	M_offsetCol + this->M_mesh.getFractureFacetsIdsVector()[juncture_it.second[counter]].getIdAsCell(),
                    -QFf);  // Triplet
            } // for
        } // for
    } // for
} // StiffMatrix::assembleFracture

void StiffMatrix::assembleFractureBC()
{
    for (auto& edge_it : this->M_mesh.getBorderTipEdgesIdsVector())
    {
    	bool isD = false;
        UInt borderId = 0;

        // select which BC to apply
        for(auto border_it : edge_it.getBorderIds())
        {
        	// BC = D > N && the one with greatest id
        	if(M_bc.getBordersBCMap().at(border_it).getBCType() == Dirichlet)
        	{
        		if(!isD)
        		{
        			isD = true;
        			borderId = border_it;
        		}
        		else
        		{
        			borderId = (border_it > borderId) ? border_it : borderId;
        		}
        	}
        	else if(!isD && M_bc.getBordersBCMap().at(border_it).getBCType() == Neumann)
        	{
        		borderId = (border_it > borderId) ? border_it : borderId;
        	}
        }

        // loop over the fracture facets
        for(auto facet_it : edge_it.getSeparatedFacetsIds())
        {
        	if(this->M_mesh.getFacetsVector()[facet_it].isFracture())
        	{
        		UInt neighborIdAsCell = M_mesh.getFacetsVector()[facet_it].getFractureFacetId() + M_mesh.getCellsVector().size();
        		Real aperture = M_properties.getProperties(this->M_mesh.getFacetsVector()[facet_it].getZoneCode()).M_aperture;

                if(M_bc.getBordersBCMap().at(borderId).getBCType() == Neumann)
                {
                    // Nella cond di bordo c'è già il contributo della permeabilità e mobilità (ma non la densità!)
                    const Real Q1o = M_bc.getBordersBCMap().at(borderId).getBC()(edge_it.getCentroid()) * edge_it.getSize() * aperture;

                    M_b->operator()(M_offsetRow + neighborIdAsCell) += Q1o;
                } // if
                else if(M_bc.getBordersBCMap().at(borderId).getBCType() == Dirichlet)
                {
                    const Real alpha1 = findAlpha (facet_it, &edge_it);
                    const Real alpha2 = findDirichletAlpha (facet_it, &edge_it);

                    const Real T12 = alpha1*alpha2/(alpha1 + alpha2);
                    const Real Q12 = T12 * M_properties.getMobility();
                    const Real Q1o = Q12 * M_bc.getBordersBCMap().at(borderId).getBC()(edge_it.getCentroid());

                    M_matrixElements.emplace_back(M_offsetRow + neighborIdAsCell, M_offsetCol + neighborIdAsCell, Q12); // Triplet
                    M_b->operator()(M_offsetRow + neighborIdAsCell) += Q1o;
                } // else if
        	} // if
        }// for
    } // for
} // StiffMatrix::assembleFractureBC

} // namespace FVCode3D
