/*!
 * @file Stiffness.cpp
 * @brief This class build a Stiffness-matrix of the Darcy problem (definitions).
 */

#include <FVCode3D/assembler/Stiffness.hpp>
#include <FVCode3D/property/Properties.hpp>
#include <FVCode3D/property/Permeability.hpp>

#include <FVCode3D/export/ExportVTU.hpp>
#include <Eigen/LU>
#include <unsupported/Eigen/SparseExtra>
// MFD: Lumped version of the stiffness matrix
//#define DIAGONALZ

// MFD: compute the exact inverse of M
#define INVERSEM
// MFD: compute M^-1 C with tensor-vector product: INVERSEM needed!
#define TVP

// MFD: export some matrices, requires both INVERSEM and APPROXM
//#define MFD_VERBOSE

// Assure that INVERSEM is defined
#ifdef TVP
#define INVERSEM
#endif // TVP

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
    const Point3D borderCenter = getBorderCenter(fj);
    const Point3D cellCenter = this->M_mesh.getFractureFacetsIdsVector()[n_Id].getCentroid();

    const Real A = M_properties.getProperties(this->M_mesh.getFractureFacetsIdsVector()[n_Id].getZoneCode()).M_aperture *
                  (this->M_mesh.getNodesVector()[fj.second]-this->M_mesh.getNodesVector()[fj.first]).norm();

    const PermPtr_Type & k =
            M_properties.getProperties(this->M_mesh.getFractureFacetsIdsVector()[n_Id].getZoneCode()).M_permeability;

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

Real StiffMatrix::findAlpha (const UInt & cellId, const Facet_ID * facet) const
{
    const Point3D facetCenter = facet->getCentroid();
    const Point3D cellCenter = this->M_mesh.getCellsVector()[cellId].getCentroid();

    const PermPtr_Type & k =
            M_properties.getProperties(this->M_mesh.getCellsVector()[cellId].getZoneCode()).M_permeability;

    Point3D f = cellCenter - facetCenter;
    const Real D = sqrt(dotProduct(f, f));
    assert( D != 0 );

    f /= D;

    const Point3D normal = facet->getUnsignedNormal();

    const Real KnDotF = fabs(dotProduct(k * normal, f)); // K n * f

    return facet->getSize() * KnDotF / D;
}

Real StiffMatrix::findAlpha (const UInt & facetId, const Edge_ID * edge) const
{
    const Point3D edgeCenter = edge->getCentroid();
    const Point3D facetCenter = this->M_mesh.getFacetsVector()[facetId].getCentroid();

    const Real A = edge->getSize() *
                   M_properties.getProperties(this->M_mesh.getFacetsVector()[facetId].getZoneCode()).M_aperture;

    const PermPtr_Type & k =
            M_properties.getProperties(this->M_mesh.getFacetsVector()[facetId].getZoneCode()).M_permeability;

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

Real StiffMatrix::findDirichletAlpha (const UInt & cellId, const Facet_ID * facet) const
{
    const Point3D facetCenter = facet->getCentroid();
    const Point3D cellCenter = this->M_mesh.getCellsVector()[cellId].getCentroid();

    const PermPtr_Type & k =
            M_properties.getProperties(this->M_mesh.getCellsVector()[cellId].getZoneCode()).M_permeability;

    Point3D f = cellCenter - facetCenter;
    const Real D = sqrt(dotProduct(f, f));
    assert( D != 0 );

    f /= D;

    const Point3D normal = facet->getUnsignedNormal();

    const Real KnDotF = fabs(dotProduct(k * normal, f)); // K n * f

    return facet->getSize() * KnDotF / D;
}

Real StiffMatrix::findDirichletAlpha (const UInt & facetId, const Edge_ID * edge) const
{
    const Point3D borderCenter = edge->getCentroid();
    const Point3D facetCenter = this->M_mesh.getFacetsVector()[facetId].getCentroid();

    const Real A = edge->getSize() *
                   M_properties.getProperties(this->M_mesh.getFacetsVector()[facetId].getZoneCode()).M_aperture;

    const PermPtr_Type & k =
            M_properties.getProperties(this->M_mesh.getFacetsVector()[facetId].getZoneCode()).M_permeability;

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

    if( !this->M_mesh.getBorderFacetsIdsVector().empty() )
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
        assert( ( alpha1 + alpha2 ) != 0 );

        const Real T12 = alpha1*alpha2/(alpha1 + alpha2);
        const Real Q12 = T12 * M_properties.getMobility();

        M_matrixElements.emplace_back(M_offsetRow + neighbor1id, M_offsetCol + neighbor1id, Q12); // Triplet
        M_matrixElements.emplace_back(M_offsetRow + neighbor1id, M_offsetCol + neighbor2id, -Q12); // Triplet
        M_matrixElements.emplace_back(M_offsetRow + neighbor2id, M_offsetCol + neighbor1id, -Q12); // Triplet
        M_matrixElements.emplace_back(M_offsetRow + neighbor2id, M_offsetCol + neighbor2id, Q12); // Triplet
    } // for

    for (auto& facet_it : this->M_mesh.getFractureFacetsIdsVector())
    {
        const PermPtr_Type & F_permeability = M_properties.getProperties(facet_it.getZoneCode()).M_permeability;
        const Real F_aperture = M_properties.getProperties(facet_it.getZoneCode()).M_aperture;
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
        const Real Q1f = T1f * M_properties.getMobility();

        assert( alphaf != 0 && alpha2 != 0 );
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

            M_b->operator()(M_offsetRow + neighbor1id) -= Q1o;
        } // if
        else if(M_bc.getBordersBCMap().at(borderId).getBCType() == Dirichlet)
        {
            const Real alpha1 = findAlpha (neighbor1id, &facet_it);
            const Real alpha2 = findDirichletAlpha (neighbor1id, &facet_it);

            assert( alpha1 != 0 && alpha2 != 0 );
            const Real T12 = 2 * alpha1*alpha2/(alpha1 + alpha2);
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

            assert( a_sum != 0 );

            for(UInt counter = 0; counter < alphas.size(); ++counter)
            {
                const Real QFf = alphaF * alphas[counter] * M_properties.getMobility() / a_sum;
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
                const UInt neighborIdAsCell = M_mesh.getFacetsVector()[facet_it].getFractureFacetId() + M_mesh.getCellsVector().size();
                const Real aperture = M_properties.getProperties(this->M_mesh.getFacetsVector()[facet_it].getZoneCode()).M_aperture;

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
                    const Real Q12 = T12 * M_properties.getMobility();
                    const Real Q1o = Q12 * M_bc.getBordersBCMap().at(borderId).getBC()(edge_it.getCentroid());

                    M_matrixElements.emplace_back(M_offsetRow + neighborIdAsCell, M_offsetCol + neighborIdAsCell, Q12); // Triplet
                    M_b->operator()(M_offsetRow + neighborIdAsCell) += Q1o;
                } // else if
            } // if
        }// for
    } // for
} // StiffMatrix::assembleFractureBC

void StiffMatrix::assembleMFD()
{
    std::cout<<std::endl;
    using Eigen::Dynamic;
    using Eigen::RowMajor;
    using Eigen::ColMajor;

    std::vector<Triplet> Matrix_elements;
    std::vector<UInt> ZMatrix_elements;
    std::vector<UInt> BMatrix_elements;
    std::vector<UInt> CMatrix_elements;
    const Real tCoeff=2.;

    // Local Matrices for Mimetic
    Eigen::Matrix<Real,Dynamic,3> Np;         // facet normals
    Eigen::Matrix<Real,3,3> Kp;               // permability tensor
    Eigen::Matrix<Real,1,Dynamic> Bp;         // areas
    Eigen::Matrix<Real,1,Dynamic> Cp;         // areas no neumann/fractures
    Eigen::Matrix<Real,Dynamic,3> Rp;         // centroid * area
    Eigen::Matrix<Real,Dynamic,Dynamic> Z0p;  // Component of Z matrix (similar to Mp0)
    Eigen::Matrix<Real,Dynamic,Dynamic> Z1p;  // Component of Z Matrix (similar to Mp1)
    Eigen::Matrix<Real,Dynamic,Dynamic> Zp;   // Z inverse matrix for local internal product= \f$ Mp^{-1}\f$
    Eigen::Matrix<Real,Dynamic,Dynamic> M0p;  // Component of M matrix (consistency)
    Eigen::Matrix<Real,Dynamic,Dynamic> M1p;  // Component of M Matrix (stability)
    Eigen::Matrix<Real,Dynamic,Dynamic> Mp;   // M matrix for internal product
    Eigen::Matrix<Real,Dynamic,Dynamic> Sp;   // Base for the column space of Np
    Eigen::Matrix<Real,Dynamic,Dynamic> Qp;   // Base for the column space of Rp

    std::vector<UInt> globalNeumannFaces;     // Indices of Neumann/fracture facets

    // Extract info of mesh. auto&& resolves const &
    auto& cellVectorRef  = this->M_mesh.getCellsVector();
    auto& facetVectorRef = this->M_mesh.getFacetsVector();
    auto& fracVectorRef = this->M_mesh.getFractureFacetsIdsVector();
    const UInt numFacets = facetVectorRef.size();
    const UInt numCells = cellVectorRef.size();
    const UInt numFracs = fracVectorRef.size();
    const UInt numCellsTot = numCells + numFracs;

    // Sizing global matrices
    Eigen::SparseMatrix<Real, RowMajor> Z;  // Z inverse matrix for internal product= \f$ M^{-1}\f$
    Eigen::SparseMatrix<Real, RowMajor> M;  // M matrix for internal product
#ifdef TVP
    Eigen::SparseMatrix<Real, ColMajor> Bt;  // Bt matrix for signed area
#else // TVP
    Eigen::SparseMatrix<Real, RowMajor> B;  // B matrix for signed area
#endif // TVP
    Eigen::SparseMatrix<Real, RowMajor> C;  // C matrix for signed area, no neumann/fractures
    Eigen::SparseMatrix<Real, RowMajor> Bd; // B matrix for Dirichlet area
    Eigen::SparseMatrix<Real, ColMajor> T(numCells,numCells);
    Eigen::SparseMatrix<Real, ColMajor> Td(numCells,numFacets); // T matrix that contains the Dirichlet contributes
    ZMatrix_elements.resize(numFacets,0);
#ifdef TVP
    BMatrix_elements.resize(numFacets,0);
#else // TVP
    BMatrix_elements.resize(numCells,0);
#endif // TVP
    CMatrix_elements.resize(numCells,0);
    Z.resize(numFacets,numFacets);
    M.resize(numFacets,numFacets);
#ifdef TVP
    Bt.resize(numFacets,numCells);
#else // TVP
    B.resize(numCells,numFacets);
#endif // TVP
    C.resize(numCells,numFacets);
    Bd.resize(numFacets,numFacets);

    // First loop to size matrices and avoid memory realloc
    for (auto& cell : cellVectorRef)
    {
        const std::vector<UInt> & cellFacetsId( cell.getFacetsIds() );
        UInt numCellFacets = cellFacetsId.size();
        UInt cellId        = cell.getId();

        for(UInt localFacetId=0; localFacetId<numCellFacets; ++localFacetId)
        {
            const UInt globalFacetId = cellFacetsId[localFacetId];
            const Rigid_Mesh::Facet & fac = facetVectorRef[globalFacetId];
#ifdef TVP
            BMatrix_elements[globalFacetId] += 1;
#else // TVP
            BMatrix_elements[cellId] += 1;
#endif // TVP

            if(
                ( fac.isBorderFacet() &&
                  (M_bc.getBordersBCMap().at(fac.getBorderId()).getBCType() == Neumann)
                ) ||
                ( fac.isFracture() )
              )
            {
                continue;
            }

            CMatrix_elements[cellId] += 1;
#ifdef DIAGONALZ
            ZMatrix_elements[globalFacetId] += 1;
#else // DIAGONALZ
            ZMatrix_elements[globalFacetId] += numCellFacets;
#endif // DIAGONALZ
        }
    }
    Z.reserve(ZMatrix_elements);
    M.reserve(ZMatrix_elements);
#ifdef TVP
    Bt.reserve(BMatrix_elements);
#else // TVP
    B.reserve(BMatrix_elements);
#endif // TVP
    C.reserve(CMatrix_elements);
    ZMatrix_elements.clear();
    ZMatrix_elements.shrink_to_fit();
    BMatrix_elements.clear();
    BMatrix_elements.shrink_to_fit();
    CMatrix_elements.clear();
    CMatrix_elements.shrink_to_fit();

    // Loop on cells
    for (auto& cell : cellVectorRef)
    {
        const std::vector<UInt> & cellFacetsId( cell.getFacetsIds() );
        const UInt numCellFacets  = cellFacetsId.size();
        auto& K = M_properties.getProperties(cell.getZoneCode()).M_permeability;
        auto cellBaricenter = cell.getCentroid();
        auto cellMeasure    = cell.getVolume();
        auto cellId         = cell.getId();

        if(cellId % 500 == 0)
        {
            std::cout<<"Done "<< cellId<<" Cells"<<std::endl;
        }

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
        Cp.resize(Eigen::NoChange,numCellFacets);
        Cp.setZero();

        Kp(0,0) = K->operator()(0,0);
        Kp(0,1) = K->operator()(0,1);
        Kp(0,2) = K->operator()(0,2);
        Kp(1,0) = K->operator()(1,0);
        Kp(1,1) = K->operator()(1,1);
        Kp(1,2) = K->operator()(1,2);
        Kp(2,0) = K->operator()(2,0);
        Kp(2,1) = K->operator()(2,1);
        Kp(2,2) = K->operator()(2,2);

        for(UInt localFacetId=0; localFacetId<numCellFacets; ++localFacetId)
        {
            const UInt globalFacetId = cellFacetsId[localFacetId];
            const Rigid_Mesh::Facet & fac = facetVectorRef[globalFacetId];

            const Point3D facetBaricenter = fac.getCentroid();
            const Point3D facetNormal     = fac.getUnsignedNormal();
            Point3D g   = facetBaricenter - cellBaricenter;
            const Real facetMeasure       = fac.area();
            const Real dotp   = dotProduct(g,facetNormal);
            const Real alpha  = (dotp >=0. ? 1.0 : -1.0);

            //! @todo I use the formulation of Nicola (to be reviewed)
            // BEWARE FOR THE CONDITION ON VELOCITY I NEED THE AVERAGE
            // VELOCITY NOT THE FLUX
            Np(localFacetId,0)    = facetNormal[0];
            Np(localFacetId,1)    = facetNormal[1];
            Np(localFacetId,2)    = facetNormal[2];
            Np *= Kp;

            Bp(0,localFacetId)    = alpha * facetMeasure;

            g[0] = std::fabs(g[0]) < 1e-15 ? 0 : g[0]; // TODO set relative eps
            g[1] = std::fabs(g[1]) < 1e-15 ? 0 : g[1]; // TODO set relative eps
            g[2] = std::fabs(g[2]) < 1e-15 ? 0 : g[2]; // TODO set relative eps

            Rp(localFacetId,0) = alpha * g[0] * facetMeasure;
            Rp(localFacetId,1) = alpha * g[1] * facetMeasure;
            Rp(localFacetId,2) = alpha * g[2] * facetMeasure;

            if(
                ( fac.isBorderFacet() &&
                  (M_bc.getBordersBCMap().at(fac.getBorderId()).getBCType() == Neumann)
                ) ||
                ( fac.isFracture() )
              )
            {
                globalNeumannFaces.push_back(globalFacetId);
                continue;
            }

            Cp(0,localFacetId)    = alpha * facetMeasure;

            if (fac.isBorderFacet() && (M_bc.getBordersBCMap().at(fac.getBorderId()).getBCType() == Dirichlet))
            {
                Bd.coeffRef(globalFacetId,globalFacetId) = alpha * facetMeasure;
            }
        }

#ifdef INVERSEM
        Sp = Np.fullPivLu().image(Np);

        Eigen::Matrix<Real,Dynamic,Dynamic> NtN = Sp.transpose() * Sp;
        Eigen::Matrix<Real,Dynamic,Dynamic> NNtNiNt = Sp * NtN.inverse() * Sp.transpose();

        M0p.resize(numCellFacets,numCellFacets);
        M0p.setZero();
        M1p.resize(numCellFacets,numCellFacets);
        M1p.setZero();
        Mp.resize(numCellFacets,numCellFacets);
        Mp.setZero();

        M0p = ( 1. / cellMeasure ) *
              (Rp * Kp.inverse() * Rp.transpose());
        M1p = M0p.trace() * tCoeff / numCellFacets *
              ( Eigen::MatrixXd::Identity(numCellFacets,numCellFacets) -
                NNtNiNt
              );

        Mp  = M0p + M1p;
#else // INVERSEM
        Qp = Rp.fullPivLu().image(Rp);

        Eigen::Matrix<Real,Dynamic,Dynamic> RtR = Qp.transpose() * Qp;
        Eigen::Matrix<Real,Dynamic,Dynamic> RRtRiRt = Qp * RtR.inverse() * Qp.transpose();

        Z0p = ( 1. / cellMeasure ) *
              (Np * Kp.inverse() * Np.transpose());
        Z1p = Z0p.trace() * tCoeff / numCellFacets *
              ( Eigen::MatrixXd::Identity(numCellFacets,numCellFacets) -
                RRtRiRt
//                Qp * Qp.transpose()
              );

        Zp  = Z0p + Z1p;

#ifdef MFD_VERBOSE
        std::stringstream  ss;
        ss << "./mMatrix/RRtRiRt_" << cellId << ".m";
        Eigen::Matrix<Real,3,3> RtR = Rp.transpose() * Rp;
        Real det = RtR.determinant();
        std::cout<<"ID: "<<cellId<<" Det: "<<det<<std::endl;
        if(std::fabs(det) >= 1e-10)
        {
            Eigen::Matrix<Real,Dynamic,Dynamic> RRtRiRt = Rp * RtR.inverse() * Rp.transpose();
            Eigen::saveMarket( RRtRiRt, ss.str() );
        }
        Eigen::Matrix<Real,Dynamic,Dynamic> QQt = Qp * Qp.transpose();
        ss.str(std::string());
        ss << "./mMatrix/QQt_" << cellId << ".m";
        Eigen::saveMarket( QQt, ss.str() );
        ss.str(std::string());
        ss << "./mMatrix/Rp_" << cellId << ".m";
        Eigen::saveMarket( Rp, ss.str() );
        ss.str(std::string());
        ss << "./mMatrix/Qp_" << cellId << ".m";
        Eigen::saveMarket( Qp, ss.str() );
#endif // MFD_VERBOSE
#endif // INVERSEM

        for(UInt iloc=0; iloc<numCellFacets; ++iloc)
        {
            const UInt i = cellFacetsId[iloc];
            const Real Bcoeff = Bp(0,iloc);
            const Real Ccoeff = Cp(0,iloc);

            if(Bcoeff != 0.)
            {
#ifdef TVP
                Bt.insert(i,cellId) = Bcoeff;
#else // TVP
                B.insert(cellId,i) = Bcoeff;
#endif // TVP
                if(Ccoeff != 0.)
                {
                    C.insert(cellId,i) = Ccoeff;
                }
            }

#ifdef INVERSEM
            Real Mcoeff = Mp(iloc,iloc);

            if( Mcoeff != 0. )
            {
                M.coeffRef(i,i) += Mcoeff;
            }

            for(UInt jloc=iloc+1; jloc<numCellFacets; ++jloc)
            {
                const UInt j = cellFacetsId[jloc];
                Mcoeff = Mp(iloc,jloc);

                if (Mcoeff != 0.)
                {
#ifdef DIAGONALZ
                    M.coeffRef(i,i) += Mcoeff;
                    M.coeffRef(j,j) += Mcoeff;
#else // DIAGONALZ
                    M.coeffRef(i,j) += Mcoeff;
                    M.coeffRef(j,i) += Mcoeff;
#endif // DIAGONALZ
                }
            }
#else //INVERSEM
            Real Zcoeff = Zp(iloc,iloc);

            if( Zcoeff != 0. )
            {
                Z.coeffRef(i,i) += Zcoeff;
            }

            for(UInt jloc=iloc+1; jloc<numCellFacets; ++jloc)
            {
                const UInt j = cellFacetsId[jloc];
                Zcoeff = Zp(iloc,jloc);

                if (Zcoeff != 0.)
                {
#ifdef DIAGONALZ
                    Z.coeffRef(i,i) += Zcoeff;
                    Z.coeffRef(j,j) += Zcoeff;
#else // DIAGONALZ
                    Z.coeffRef(i,j) += Zcoeff;
                    Z.coeffRef(j,i) += Zcoeff;
#endif // DIAGONALZ
                }
            }
#endif // INVERSEM
        }
    }

#ifdef INVERSEM
#ifdef TVP
    Eigen::SparseMatrix<Real, ColMajor> WtCM;
    WtCM.resize(numFacets,numCells);
    std::cout << "Compute B M^-1" << std::endl;

    Eigen::ConjugateGradient< Eigen::SparseMatrix<Real, RowMajor> > cg;

    cg.setMaxIterations(10e3);
    cg.setTolerance(1e-6);

    cg.compute(M);
    for(UInt i=0; i<Bt.cols(); ++i)
    {
        if(i % 100 == 0)
        {
            std::cout << "System: " << i << " ..." << std::endl;
        }

        const Vector t_rhs(Bt.col(i));
        Vector t_x;
        t_x = cg.solve(t_rhs);
        WtCM.col(i) = t_x.sparseView();

#ifdef MFD_VERBOSE
        if(i % 100 == 0)
        {
            std::cout << "# It:  " << cg.iterations() << std::endl;
            std::cout << "# Err: " << cg.error() << std::endl;
        }
#endif // MFD_VERBOSE
    }
    std::cout << "System: " << (Bt.cols() - 1) << " done." << std::endl;
    Eigen::SparseMatrix<Real, RowMajor> Wt(WtCM);
    WtCM.resize(0,0);
    WtCM.data().squeeze();
    Eigen::SparseMatrix<Real, RowMajor> BtRM(Bt);
    Bt.resize(0,0);
    Bt.data().squeeze();

    std::cout << "Zero the Neumann row of W" << std::endl;
    for (UInt i = 0; i < globalNeumannFaces.size(); ++i)
    {
        for (Eigen::SparseMatrix<Real, RowMajor>::InnerIterator it(Wt,globalNeumannFaces[i]); it; ++it)
        {
            it.valueRef() = 0.;
        }
        for (Eigen::SparseMatrix<Real, RowMajor>::InnerIterator it(BtRM,globalNeumannFaces[i]); it; ++it)
        {
            Wt.coeffRef(globalNeumannFaces[i], it.col()-1) = it.valueRef();
        }
    }
    Wt.prune(0.);
    std::cout<<" Matrix Wt has "<<Wt.rows()<<" rows and "<<Wt.cols()<<" columns and "<<Wt.nonZeros()<<" Non zeros"<<std::endl;
#else // TVP
    Eigen::Matrix<Real,Dynamic,Dynamic> invM(M);
    std::cout<<"Compute inverse of M"<<std::endl;
    invM = invM.inverse();

    std::cout<<"Zero the Neumann row of M"<<std::endl;
    for(UInt i = 0; i < globalNeumannFaces.size(); ++i)
    {
        invM.row(globalNeumannFaces[i]).setZero();
        invM.coeffRef(globalNeumannFaces[i], globalNeumannFaces[i]) = 1.;
    }
    Eigen::SparseMatrix<Real, RowMajor> invMSp;  // M matrix for internal product
    invMSp = invM.sparseView();
#endif //TVP
#else // INVERSEM
    for (UInt i = 0; i < globalNeumannFaces.size(); ++i)
    {
        for (Eigen::SparseMatrix<Real, RowMajor>::InnerIterator it(Z,globalNeumannFaces[i]); it; ++it)
        {
            it.valueRef() = 0.;
        }
        Z.coeffRef(globalNeumannFaces[i], globalNeumannFaces[i]) = 1.;
    }

//    for(UInt i = 0; i < globalNeumannFaces.size(); ++i)
//    {
//        for(UInt col = 0; col < Z.cols(); ++col)
//        {
//            Z.coeffRef(globalNeumannFaces[i],col) = 0.;
//        }
//        Z.coeffRef(globalNeumannFaces[i], globalNeumannFaces[i]) = 1.;
//    }
#endif // INVERSEM

    // Now I Need to build T
    std::cout<<"Building Trasmissibility matrix (global)"<<std::endl;

#ifdef INVERSEM
#ifdef TVP
    T = Wt.transpose() * C.transpose();
#else // TVP
    T = B * invMSp * C.transpose();
#endif // TVP
#else // INVERSEM
    T = B * Z * C.transpose();
#endif // INVERSEM
    T.makeCompressed();

    std::cout<<"Building Trasmissibility Dirichelt matrix (global)"<<std::endl;

#ifdef INVERSEM
#ifdef TVP
    Td = Wt.transpose() * Bd.transpose();
#else // TVP
    Td = B * invMSp * Bd.transpose();
#endif // TVP
#else // INVERSEM
    Td = B * Z * Bd.transpose();
#endif // INVERSEM
    Td.makeCompressed();

    // Clean up to seve memory (I hope)
    std::cout<<" Matrix T has "<<T.rows()<<" rows and "<<T.cols()<<" columns and "<<T.nonZeros()<<" Non zeros"<<std::endl;
    std::cout<<" Matrix Td has "<<Td.rows()<<" rows and "<<Td.cols()<<" columns and "<<Td.nonZeros()<<" Non zeros"<<std::endl;
    std::cout<<" Num Cells: "<<numCells<<std::endl;
    std::cout<<" Num Facets: "<<numFacets<<std::endl;
    std::cout<<" Num Cells + Fracs: "<<numCellsTot<<std::endl;
#ifdef MFD_VERBOSE
    std::cout.flush();
    Eigen::saveMarket( T, "./mMatrix/T.m" );
    Eigen::saveMarket( Td, "./mMatrix/Td.m" );
#ifdef TVP
    Eigen::saveMarket( BtRM, "./mMatrix/Bt.m" );
#else // TVP
    Eigen::saveMarket( B, "./mMatrix/B.m" );
#endif // TVP
    Eigen::saveMarket( C, "./mMatrix/C.m" );
    Eigen::saveMarket( Bd, "./mMatrix/Bd.m" );
#ifdef INVERSEM
    Eigen::saveMarket( invM, "./mMatrix/invM.m" );
#else // INVERSEM
    Eigen::saveMarket( Z, "./mMatrix/Z.m" );
#endif // INVERSEM
#endif // MFD_VERBOSE

    // assemble fractures with FV
    for (auto& facet_it : this->M_mesh.getFractureFacetsIdsVector())
    {
        auto& F_permeability = M_properties.getProperties(facet_it.getZoneCode()).M_permeability;
        const Real F_aperture = M_properties.getProperties(facet_it.getZoneCode()).M_aperture;
        const Real Df = F_aperture/2.;
        const Point3D normal = facet_it.getUnsignedNormal();
        const Real KnDotF = fabs(dotProduct(F_permeability * normal, normal));
        const Real alphaf = facet_it.getSize() * KnDotF / Df;

        const UInt neighbor1id = facet_it.getSeparatedCellsIds()[0];
        const UInt neighbor2id = facet_it.getSeparatedCellsIds()[1];

        const Real alpha1 = findAlpha(neighbor1id, &facet_it);
        const Real alpha2 = findAlpha(neighbor2id, &facet_it);

        const Real T1f = alpha1*alphaf/(alpha1 + alphaf);
        const Real Q1f = T1f * M_properties.getMobility();

        const Real T2f = alphaf*alpha2/(alphaf + alpha2);
        const Real Q2f = T2f * M_properties.getMobility();

        Matrix_elements.emplace_back(Triplet (neighbor1id, neighbor1id, Q1f));
        Matrix_elements.emplace_back(Triplet (neighbor1id, facet_it.getIdAsCell(), -Q1f));
        Matrix_elements.emplace_back(Triplet (facet_it.getIdAsCell(), neighbor1id, -Q1f));
        Matrix_elements.emplace_back(Triplet (facet_it.getIdAsCell(), facet_it.getIdAsCell(), Q1f));

        Matrix_elements.emplace_back(Triplet (neighbor2id, neighbor2id, Q2f));
        Matrix_elements.emplace_back(Triplet (neighbor2id, facet_it.getIdAsCell(), -Q2f));
        Matrix_elements.emplace_back(Triplet (facet_it.getIdAsCell(), neighbor2id, -Q2f));
        Matrix_elements.emplace_back(Triplet (facet_it.getIdAsCell(), facet_it.getIdAsCell(), Q2f));

        for (auto juncture_it : facet_it.getFractureNeighbors())
        {
            std::vector<Real> alphas;
            const Real alphaF = findFracturesAlpha (juncture_it.first, facet_it.getFractureId());

            for (auto neighbors_it : juncture_it.second)
            {
                alphas.emplace_back(findFracturesAlpha (juncture_it.first, neighbors_it));
            }

            Real a_sum = alphaF;
            auto sum_maker = [&a_sum](Real item){a_sum += item;};
            std::for_each(alphas.begin(), alphas.end(), sum_maker);

            for (UInt counter = 0; counter < alphas.size(); ++counter)
            {
                const Real QFf = alphaF * alphas[counter] * M_properties.getMobility() / a_sum;
                Matrix_elements.emplace_back(Triplet (facet_it.getIdAsCell(), facet_it.getIdAsCell(), QFf));
                Matrix_elements.emplace_back(Triplet (facet_it.getIdAsCell(),
                this->M_mesh.getFractureFacetsIdsVector()[juncture_it.second[counter]].getIdAsCell(), -QFf));
            }
            alphas.clear();
        }
    }

    for(UInt i=0; i<this->M_size; ++i)
    {
        M_b->operator()(i) = 0.;
    }

    // assemble BC on porous matrix
    for (auto& facet_it : this->M_mesh.getBorderFacetsIdsVector())
    {
        const UInt neighbor1id = facet_it.getSeparatedCellsIds()[0];
        const UInt borderId = facet_it.getBorderId();
        const UInt facetId = facet_it.getId();

        if(M_bc.getBordersBCMap().at(borderId).getBCType() == Neumann )
        {
            // Nella cond di bordo c'è già il contributo della permeabilità e mobilità (ma non la densità!)
            const Real Q1o = M_bc.getBordersBCMap().at(borderId).getBC()(facet_it.getCentroid()) * facet_it.getSize();

            M_b->operator()(neighbor1id) -= Q1o;
        }
        else if(M_bc.getBordersBCMap().at(borderId).getBCType() == Dirichlet)
        {
            for (Eigen::SparseMatrix<Real>::InnerIterator it(Td,facetId); it; ++it)
            {
                M_b->operator()(it.row()) += it.value() * M_bc.getBordersBCMap().at(borderId).getBC()(facet_it.getCentroid());
            }
        }
    }

    // assemble BC on fractures
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
                const Real aperture = M_properties.getProperties(this->M_mesh.getFacetsVector()[facet_it].getZoneCode()).M_aperture;

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

                    const Real T12 = 2 * alpha1*alpha2/(alpha1 + alpha2);
                    const Real Q12 = T12 * M_properties.getMobility();
                    const Real Q1o = Q12 * M_bc.getBordersBCMap().at(borderId).getBC()(edge_it.getCentroid());

                    Matrix_elements.emplace_back(M_offsetRow + neighborIdAsCell, M_offsetCol + neighborIdAsCell, Q12); // Triplet
                    M_b->operator()(M_offsetRow + neighborIdAsCell) += Q1o;
                } // else if
            } // if
        }// for
    } // for

    // Put everything together
    Eigen::SparseMatrix<Real,ColMajor> Mborder(numCellsTot,numCellsTot);
    Mborder.setFromTriplets(Matrix_elements.begin(), Matrix_elements.end());
    this->M_Matrix->resize(numCellsTot,numCellsTot);
    T.conservativeResize(numCellsTot,numCellsTot);
    *(this->M_Matrix) = T + Mborder;
//    Eigen::saveMarket( *(this->M_Matrix), "./mMatrix/A.m" );
//    Eigen::saveMarket( *M_b, "./mMatrix/RHS.m" );
    std::cout<<" Assembling ended"<<std::endl;
} // StiffMatrix::assembleMFD

} // namespace FVCode3D
