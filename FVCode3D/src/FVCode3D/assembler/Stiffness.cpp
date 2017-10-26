/*!
 * @file Stiffness.cpp
 * @brief This class build a Stiffness-matrix of the Darcy problem (definitions).
 */

#include <FVCode3D/assembler/Stiffness.hpp>
#include <FVCode3D/property/Properties.hpp>
#include <FVCode3D/property/Permeability.hpp>
#include <FVCode3D/assembler/global_operator.hpp>

#include <FVCode3D/export/ExportVTU.hpp>
#include <Eigen/LU>
#include <unsupported/Eigen/SparseExtra>

// MFD: export some matrices, requires both INVERSEM and APPROXM
//#define MFD_VERBOSE

namespace FVCode3D
{

Point3D StiffMatrixFV::getBorderCenter(Fracture_Juncture fj) const
{
    return (this->M_mesh.getNodesVector()[fj.first] +
                (this->M_mesh.getNodesVector()[fj.second] -
                 this->M_mesh.getNodesVector()[fj.first]
                )/2.
           );
}

Real StiffMatrixFV::findFracturesAlpha (const Fracture_Juncture fj, const UInt n_Id) const
{
    const Point3D borderCenter = getBorderCenter(fj);
    const Point3D cellCenter = this->M_mesh.getFractureFacetsIdsVector()[n_Id].getCentroid();

    const Real A = M_mesh.getPropertiesMap().getProperties(this->M_mesh.getFractureFacetsIdsVector()[n_Id].getZoneCode()).M_aperture *
                  (this->M_mesh.getNodesVector()[fj.second]-this->M_mesh.getNodesVector()[fj.first]).norm();

    const PermPtr_Type & k =
            M_mesh.getPropertiesMap().getProperties(this->M_mesh.getFractureFacetsIdsVector()[n_Id].getZoneCode()).M_permeability;

    Point3D f = borderCenter - cellCenter;
    const Real D = sqrt(dotProduct(f, f));
    assert( D != 0 );

    f /= D;

    Point3D normal = crossProduct(f, this->M_mesh.getNodesVector()[fj.second]-this->M_mesh.getNodesVector()[fj.first]); // k = f x l
    normal = crossProduct(this->M_mesh.getNodesVector()[fj.second]-this->M_mesh.getNodesVector()[fj.first], normal); // n = l x k
    normal.normalize();

    const Real KnDotF = fabs(dotProduct(k * normal, f)); // K n * f

    return A * KnDotF / D;
}

Real StiffMatrixFV::findAlpha (const UInt & cellId, const Facet_ID * facet) const
{
    const Point3D facetCenter = facet->getCentroid();
    const Point3D cellCenter = this->M_mesh.getCellsVector()[cellId].getCentroid();

    const PermPtr_Type & k =
            M_mesh.getPropertiesMap().getProperties(this->M_mesh.getCellsVector()[cellId].getZoneCode()).M_permeability;

    Point3D f = cellCenter - facetCenter;
    const Real D = sqrt(dotProduct(f, f));
    assert( D != 0 );

    f /= D;

    const Point3D normal = facet->getUnsignedNormal();

    const Real KnDotF = fabs(dotProduct(k * normal, f)); // K n * f

    return facet->getSize() * KnDotF / D;
}

Real StiffMatrixFV::findAlpha (const UInt & facetId, const Edge_ID * edge) const
{
    const Point3D edgeCenter = edge->getCentroid();
    const Point3D facetCenter = this->M_mesh.getFacetsVector()[facetId].getCentroid();

    const Real A = edge->getSize() *
                   M_mesh.getPropertiesMap().getProperties(this->M_mesh.getFacetsVector()[facetId].getZoneCode()).M_aperture;

    const PermPtr_Type & k =
            M_mesh.getPropertiesMap().getProperties(this->M_mesh.getFacetsVector()[facetId].getZoneCode()).M_permeability;

    Point3D f = edgeCenter - facetCenter;
    const Real D = sqrt(dotProduct(f, f));
    assert( D != 0 );

    f /= D;

    Point3D normal = crossProduct(  f,
                            this->M_mesh.getNodesVector()[edge->getEdge().getVerticesIds()[0]] -
                            this->M_mesh.getNodesVector()[edge->getEdge().getVerticesIds()[1]] ); // k = f x l
    normal = crossProduct(  this->M_mesh.getNodesVector()[edge->getEdge().getVerticesIds()[0]] -
                            this->M_mesh.getNodesVector()[edge->getEdge().getVerticesIds()[1]],
                            normal); // n = l x k
    normal.normalize();

    const Real KnDotF = fabs(dotProduct(k * normal, f)); // K n * f

    return A * KnDotF / D;
}

Real StiffMatrixFV::findDirichletAlpha (const UInt & cellId, const Facet_ID * facet) const
{
    const Point3D facetCenter = facet->getCentroid();
    const Point3D cellCenter = this->M_mesh.getCellsVector()[cellId].getCentroid();

    const PermPtr_Type & k =
            M_mesh.getPropertiesMap().getProperties(this->M_mesh.getCellsVector()[cellId].getZoneCode()).M_permeability;

    Point3D f = cellCenter - facetCenter;
    const Real D = sqrt(dotProduct(f, f));
    assert( D != 0 );

    f /= D;

    const Point3D normal = facet->getUnsignedNormal();

    const Real KnDotF = fabs(dotProduct(k * normal, f)); // K n * f

    return facet->getSize() * KnDotF / D;
}

Real StiffMatrixFV::findDirichletAlpha (const UInt & facetId, const Edge_ID * edge) const
{
    const Point3D borderCenter = edge->getCentroid();
    const Point3D facetCenter = this->M_mesh.getFacetsVector()[facetId].getCentroid();

    const Real A = edge->getSize() *
                   M_mesh.getPropertiesMap().getProperties(this->M_mesh.getFacetsVector()[facetId].getZoneCode()).M_aperture;

    const PermPtr_Type & k =
            M_mesh.getPropertiesMap().getProperties(this->M_mesh.getFacetsVector()[facetId].getZoneCode()).M_permeability;

    Point3D f = borderCenter - facetCenter;
    const Real D = sqrt(dotProduct(f, f));
    assert( D != 0 );

    f /= D;

    Point3D normal = crossProduct(  f,
                            this->M_mesh.getNodesVector()[edge->getEdge().getVerticesIds()[0]] -
                            this->M_mesh.getNodesVector()[edge->getEdge().getVerticesIds()[1]] ); // k = f x l
    normal = crossProduct(  this->M_mesh.getNodesVector()[edge->getEdge().getVerticesIds()[0]] -
                            this->M_mesh.getNodesVector()[edge->getEdge().getVerticesIds()[1]],
                            normal); // n = l x k
    normal.normalize();

    const Real KnDotF = fabs(dotProduct(k * normal, f)); // K n * f

    return A * KnDotF / D;
}

void StiffMatrixFV::assemble()
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

    if( !this->M_mesh.getBorderFacetsIdsVector().empty() )
    {
        assemblePorousMatrixBC();
    } // if

    if( !this->M_mesh.getBorderTipEdgesIdsVector().empty() )
    {
        assembleFractureBC();
    } // if
} // assemble

void StiffMatrixFV::assemblePorousMatrix()
{
    for (auto& facet_it : this->M_mesh.getInternalFacetsIdsVector())
    {
        const UInt neighbor1id = facet_it.getSeparatedCellsIds()[0];
        const UInt neighbor2id = facet_it.getSeparatedCellsIds()[1];

        const Real alpha1 = findAlpha(neighbor1id, &facet_it);
        const Real alpha2 = findAlpha(neighbor2id, &facet_it);
        assert( ( alpha1 + alpha2 ) != 0 );

        const Real T12 = alpha1*alpha2/(alpha1 + alpha2);
        const Real Q12 = T12 * M_mesh.getPropertiesMap().getMobility();

        M_matrixElements.emplace_back(M_offsetRow + neighbor1id, M_offsetCol + neighbor1id, Q12); // Triplet
        M_matrixElements.emplace_back(M_offsetRow + neighbor1id, M_offsetCol + neighbor2id, -Q12); // Triplet
        M_matrixElements.emplace_back(M_offsetRow + neighbor2id, M_offsetCol + neighbor1id, -Q12); // Triplet
        M_matrixElements.emplace_back(M_offsetRow + neighbor2id, M_offsetCol + neighbor2id, Q12); // Triplet
    } // for

    for (auto& facet_it : this->M_mesh.getFractureFacetsIdsVector())
    {
        const PermPtr_Type & F_permeability = M_mesh.getPropertiesMap().getProperties(facet_it.getZoneCode()).M_permeability;
        const Real F_aperture = M_mesh.getPropertiesMap().getProperties(facet_it.getZoneCode()).M_aperture;
        const Real Df = F_aperture/2.;
        assert( Df != 0 );

        Point3D normal = facet_it.getUnsignedNormal();

        const Real KnDotF = fabs(dotProduct(F_permeability * normal, normal));

        const Real alphaf = facet_it.getSize() * KnDotF / Df;

        const UInt neighbor1id = facet_it.getSeparatedCellsIds()[0];
        const UInt neighbor2id = facet_it.getSeparatedCellsIds()[1];

        const Real alpha1 = findAlpha(neighbor1id, &facet_it);
        const Real alpha2 = findAlpha(neighbor2id, &facet_it);

        assert( alpha1 != 0 && alphaf != 0 );
        const Real T1f = alpha1*alphaf/(alpha1 + alphaf);
        const Real Q1f = T1f * M_mesh.getPropertiesMap().getMobility();

        assert( alphaf != 0 && alpha2 != 0 );
        const Real T2f = alphaf*alpha2/(alphaf + alpha2);
        const Real Q2f = T2f * M_mesh.getPropertiesMap().getMobility();

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

void StiffMatrixFV::assemblePorousMatrixBC()
{
    for (auto& facet_it : this->M_mesh.getBorderFacetsIdsVector())
    {
        const UInt neighbor1id = facet_it.getSeparatedCellsIds()[0];
        const UInt borderId = facet_it.getBorderId();

        if(M_bc.getBordersBCMap().at(borderId).getBCType() == Neumann )
        {
            // Nella cond di bordo c'è già il contributo della permeabilità e mobilità (ma non la densità!)
            const Real Q1o = M_bc.getBordersBCMap().at(borderId).getBC()(facet_it.getCentroid()) * facet_it.getSize();

            M_b->operator()(M_offsetRow + neighbor1id) -= Q1o;
        } // if
        else if(M_bc.getBordersBCMap().at(borderId).getBCType() == Dirichlet)
        {
            const Real alpha1 = findAlpha (neighbor1id, &facet_it);
            const Real alpha2 = findDirichletAlpha (neighbor1id, &facet_it);

            assert( alpha1 != 0 && alpha2 != 0 );
            const Real T12 = 2 * alpha1*alpha2/(alpha1 + alpha2);
            const Real Q12 = T12 * M_mesh.getPropertiesMap().getMobility();
            const Real Q1o = Q12 * M_bc.getBordersBCMap().at(borderId).getBC()(facet_it.getCentroid());

            M_matrixElements.emplace_back(M_offsetRow + neighbor1id, M_offsetCol + neighbor1id, Q12); // Triplet
            M_b->operator()(M_offsetRow + neighbor1id) += Q1o;
        } // else if
    } // for
} // StiffMatrix::assemblePorousMatrixBC

void StiffMatrixFV::assembleFracture()
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

            assert( a_sum != 0 );

            for(UInt counter = 0; counter < alphas.size(); ++counter)
            {
                const Real QFf = alphaF * alphas[counter] * M_mesh.getPropertiesMap().getMobility() / a_sum;
                M_matrixElements.emplace_back(M_offsetRow + facet_it.getIdAsCell(), M_offsetCol + facet_it.getIdAsCell(), QFf); // Triplet
                M_matrixElements.emplace_back(M_offsetRow + facet_it.getIdAsCell(),
                    M_offsetCol + this->M_mesh.getFractureFacetsIdsVector()[juncture_it.second[counter]].getIdAsCell(),
                    -QFf);  // Triplet
            } // for
        } // for
    } // for
} // StiffMatrix::assembleFracture

void StiffMatrixFV::assembleFractureBC()
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
                const UInt neighborIdAsCell = M_mesh.getFacetsVector()[facet_it].getFractureFacetId() + M_mesh.getCellsVector().size();
                const Real aperture = M_mesh.getPropertiesMap().getProperties(this->M_mesh.getFacetsVector()[facet_it].getZoneCode()).M_aperture;

                if(M_bc.getBordersBCMap().at(borderId).getBCType() == Neumann)
                {
                    // Nella cond di bordo c'è già il contributo della permeabilità e mobilità (ma non la densità!)
                    const Real Q1o = M_bc.getBordersBCMap().at(borderId).getBC()(edge_it.getCentroid()) * edge_it.getSize() * aperture;

                    M_b->operator()(M_offsetRow + neighborIdAsCell) -= Q1o;
                } // if
                else if(M_bc.getBordersBCMap().at(borderId).getBCType() == Dirichlet)
                {
                    const Real alpha1 = findAlpha (facet_it, &edge_it);
                    const Real alpha2 = findDirichletAlpha (facet_it, &edge_it);

                    assert( alpha1 != 0 && alpha2 != 0 );
                    const Real T12 = 2 * alpha1*alpha2/(alpha1 + alpha2);
                    const Real Q12 = T12 * M_mesh.getPropertiesMap().getMobility();
                    const Real Q1o = Q12 * M_bc.getBordersBCMap().at(borderId).getBC()(edge_it.getCentroid());

                    M_matrixElements.emplace_back(M_offsetRow + neighborIdAsCell, M_offsetCol + neighborIdAsCell, Q12); // Triplet
                    M_b->operator()(M_offsetRow + neighborIdAsCell) += Q1o;
                } // else if
            } // if
        }// for
    } // for
} // StiffMatrix::assembleFractureBC


void StiffMatrixMFD::assemble()
{
	// Define the degrees of fredoom
	const UInt numFracture  = M_mesh.getFractureFacetsIdsVector().size();
	const UInt numFacetsTot = M_mesh.getFacetsVector().size() + numFracture;
	const UInt numCell      = M_mesh.getCellsVector().size();
	
	// Define the global inner product
	global_InnerProduct gIP(M_mesh, dType::dFacet, dType::dFacet, numFacetsTot, numFacetsTot, M_bc);
	gIP.ShowMe();
	
	// Define the global divergence operator
	global_Div gDIV(M_mesh, dType::dCell, dType::dFacet, numCell, numFacetsTot, M_bc);
	gDIV.ShowMe();
	
	// Define the coupling conditions
	CouplingConditions coupling(M_mesh, dType::dFracture, dType::dFacet, numFracture, numFacetsTot, gIP.getMatrix());
	coupling.ShowMe();
	
	// Define the flux operator
	FluxOperator fluxOP(M_mesh, dType::dFracture, dType::dFracture, numFracture, numFracture, M_bc);
	fluxOP.ShowMe();
	
	// Define the global bulk builder
	global_BulkBuilder gBulkBuilder(M_mesh, M_bc, gIP, gDIV);
	// Reserve space for the bulk matrices
	gBulkBuilder.reserve_space(*M_Matrix);
	// Build the bulk matrices
	gBulkBuilder.build(*M_Matrix);
	// Impose bulk BCs
	gIP.ImposeBC(*M_Matrix, *M_b);
	
	// Define the fracture builder
	FractureBuilder FBuilder(M_mesh, coupling, fluxOP, gIP);
	// Reserve space for the fracture matrices
	FBuilder.reserve_space(*M_Matrix);
	// Build the fracture matrices
	FBuilder.build(*M_Matrix);
	// Impose fractures BCs
	fluxOP.ImposeBConFractures(*M_Matrix, *M_b);
}

} // namespace FVCode3D
