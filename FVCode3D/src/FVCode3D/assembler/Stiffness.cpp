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
    
    // Extract info of mesh. auto&& resolves const &
    auto& cellVectorRef    = this->M_mesh.getCellsVector();
    auto& facetVectorRef   = this->M_mesh.getFacetsVector();
    auto& fracVectorRef    = this->M_mesh.getFractureFacetsIdsVector();
    const UInt numFacets        = facetVectorRef.size();
    const UInt numCells         = cellVectorRef.size();
    const UInt numFracs         = fracVectorRef.size();
    const UInt numFacetsTot     = numFacets + numFracs;
    const UInt numCellsTot      = numCells + numFracs;
    const UInt numGDL           = numFacetsTot + numCellsTot;

    // Define global matrix and vector to reserve space
    Eigen::VectorXd & rhs                       = *(this->M_b);         // rhs  
    Eigen::SparseMatrix<Real,ColMajor> & S      = *(this->M_Matrix);    // global matrix
    std::vector<UInt> Matrix_elements(numGDL,0);                        // vector to reserve space
    
    // Problem coefficients
    const Real tCoeff     = 2.;                            // mimetic coefficient         
    const Real xsi        = 1;                             // coupling coefficient
    const Real xsi0       = ( 2. * xsi - 1. ) / 4.;        // coupling modified coefficient

    // Local MFD matrices
    Eigen::Matrix<Real,Dynamic,3> Np;         // facet normals, necessary to build Mp
    Eigen::Matrix<Real,3,3> Kp;               // permability tensor
    Eigen::Matrix<Real,1,Dynamic> Bp;         // local B
    Eigen::Matrix<Real,Dynamic,3> Rp;         // centroid * area, necessary to build Mp
    Eigen::Matrix<Real,Dynamic,Dynamic> M0p;  // Component of local M matrix (consistency)
    Eigen::Matrix<Real,Dynamic,Dynamic> M1p;  // Component of local M Matrix (stability)
    Eigen::Matrix<Real,Dynamic,Dynamic> Mp;   // local M matrix for internal product
    Eigen::Matrix<Real,Dynamic,Dynamic> Sp;   // Base for the column space of Np


    // First loop to size matrices and avoid memory realloc
    for (auto& cell : cellVectorRef){
		
        const std::vector<UInt> & cellFacetsId( cell.getFacetsIds() );
        UInt numCellFacets = cellFacetsId.size();
        UInt cellId        = cell.getId();

        for(UInt localFacetId=0; localFacetId<numCellFacets; ++localFacetId){	
			
			UInt globalFacetId = cellFacetsId[localFacetId];
			const Rigid_Mesh::Facet & fac = facetVectorRef[globalFacetId];
			
			// This is to take into account the decoupling of fractures facets
			if( fac.isFracture() ){
				
				const Point3D facetNormal = fac.getUnsignedNormal();
				const Point3D g           = fac.getCentroid() - cell.getCentroid();
				const Real dotp           = dotProduct(g,facetNormal);
				
				if( dotp < 0 )
					globalFacetId = numFacets + fac.getFractureFacetId();
				
				}

			Matrix_elements[globalFacetId] += numCellFacets + 1;        // numCellFacets for the M, the 1 for the B

            Matrix_elements[numFacetsTot + cellId] += 1;                // this is for the Bt
            
        }
    }
    
    S.reserve(Matrix_elements);

    // Loop on cells to build MFD local matrices and assemble them in MFD global matrices
    for (auto& cell : cellVectorRef){
		
        const std::vector<UInt> & cellFacetsId( cell.getFacetsIds() );
        const UInt numCellFacets  = cellFacetsId.size();
        auto& K = M_properties.getProperties(cell.getZoneCode()).M_permeability;
        const Real cellMeasure    = cell.getVolume();
        const UInt cellId         = cell.getId();
        auto cellBaricenter = cell.getCentroid();

        if(cellId % 500 == 0)
            std::cout<<"Done "<< cellId<<" Cells"<<std::endl;

        // Resize local matrices
        Np.resize(numCellFacets,Eigen::NoChange);
        Np.setZero();
        Rp.resize(numCellFacets,Eigen::NoChange);
        Rp.setZero();
        Bp.resize(Eigen::NoChange,numCellFacets);
        Bp.setZero();

        Kp(0,0) = K->operator()(0,0);
        Kp(0,1) = K->operator()(0,1);
        Kp(0,2) = K->operator()(0,2);
        Kp(1,0) = K->operator()(1,0);
        Kp(1,1) = K->operator()(1,1);
        Kp(1,2) = K->operator()(1,2);
        Kp(2,0) = K->operator()(2,0);
        Kp(2,1) = K->operator()(2,1);
        Kp(2,2) = K->operator()(2,2);

		// Loop on facets to build Rp, Np, Bp, Dtp
        for(UInt localFacetId=0; localFacetId<numCellFacets; ++localFacetId){
			
            const UInt globalFacetId = cellFacetsId[localFacetId];
            const Rigid_Mesh::Facet & fac = facetVectorRef[globalFacetId];

            const Point3D facetNormal     = fac.getUnsignedNormal();
            const Point3D facetBaricenter = fac.getCentroid();
            Point3D g                     = facetBaricenter - cellBaricenter;
            const Real facetMeasure    = fac.area();
            const Real dotp            = dotProduct(g,facetNormal);
            const Real alpha           = (dotp >=0. ? 1.0 : -1.0);

            //! @todo I use the formulation of Nicola (to be reviewed)
            // BEWARE FOR THE CONDITION ON VELOCITY I NEED THE AVERAGE
            // VELOCITY NOT THE FLUX
            Np(localFacetId,0) = facetNormal[0];
            Np(localFacetId,1) = facetNormal[1];
            Np(localFacetId,2) = facetNormal[2];
            Np *= Kp;

            Bp(0,localFacetId) =  - alpha * facetMeasure;    // The minus because I have changed the sign to the conservation equation

            g[0] = std::fabs(g[0]) < 1e-15 ? 0 : g[0]; // TODO set relative eps
            g[1] = std::fabs(g[1]) < 1e-15 ? 0 : g[1]; // TODO set relative eps
            g[2] = std::fabs(g[2]) < 1e-15 ? 0 : g[2]; // TODO set relative eps

            Rp(localFacetId,0) = alpha * g[0] * facetMeasure;
            Rp(localFacetId,1) = alpha * g[1] * facetMeasure;
            Rp(localFacetId,2) = alpha * g[2] * facetMeasure;

            if (fac.isBorderFacet() && (M_bc.getBordersBCMap().at(fac.getBorderId()).getBCType() == Dirichlet))
                rhs[globalFacetId] = - alpha * facetMeasure;

        }

		// Building the Mp matrix
        Sp = Np.fullPivLu().image(Np);    // Because otherwise NtN could be singulare in the unlucky case of 2 parallelel faces

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


        for(UInt iloc=0; iloc<numCellFacets; ++iloc){
			
            UInt i = cellFacetsId[iloc];
            const Rigid_Mesh::Facet & fac = facetVectorRef[i];
            const Real Bcoeff = Bp(0,iloc);
            
            // This is to take into account the decoupling of fractures facets
            if( fac.isFracture() ){
				
				const Point3D facetNormal = fac.getUnsignedNormal();
				const Point3D g           = fac.getCentroid() - cell.getCentroid();
				const Real dotp           = dotProduct(g,facetNormal);
				
				if( dotp < 0 )
					i = numFacets + fac.getFractureFacetId();
			}

            if( Bcoeff != 0. ){

				S.insert(numFacetsTot + cellId, i) = Bcoeff;
				S.insert(i, numFacetsTot + cellId) = Bcoeff;
			}

            Real Mcoeff = Mp(iloc,iloc);

            if( Mcoeff != 0. )
                S.coeffRef(i,i) += Mcoeff;

            for(UInt jloc=iloc+1; jloc<numCellFacets; ++jloc){

                UInt j = cellFacetsId[jloc];
                const Rigid_Mesh::Facet & fac = facetVectorRef[j];
                Mcoeff = Mp(iloc,jloc);
                
                // This is to take into account the decoupling of fractures facets
				if( fac.isFracture() ){
				
					const Point3D facetNormal = fac.getUnsignedNormal();
					const Point3D g           = fac.getCentroid() - cell.getCentroid();
					const Real dotp           = dotProduct(g,facetNormal);
				
					if( dotp < 0 )
						j = numFacets + fac.getFractureFacetId();		
				}

                if (Mcoeff != 0.){
                    S.coeffRef(i,j) += Mcoeff;
                    S.coeffRef(j,i) += Mcoeff;
                }
                
            }

        }
    }
    

    std::cout << "Zero the Neumann row of S" << std::endl;
    
	for( auto facet_it : this->M_mesh.getBorderFacetsIdsVector() ){
		
		UInt borderId = facet_it.getBorderId();
		UInt facetId  = facet_it.getId();
		
		if( this->M_bc.getBordersBCMap().at(borderId).getBCType() == Neumann ){
			
			S.prune( [facetId]( UInt i, UInt, Real ){ return ( i!= facetId ); });
			
			S.coeffRef(facetId,facetId) = 1;
			
			}
		
		}  


    std::cout<<" Num Cells: "<<numCells<<std::endl;
    std::cout<<" Num Facets: "<<numFacets<<std::endl;
    std::cout<<" Num Facets + Fracs: "<<numFacetsTot<<std::endl;
    std::cout<<" Num Cells + Fracs: "<<numCellsTot<<std::endl;
    std::cout<<" Num GDL: "<<numGDL<<std::endl;



	// Loop to avoid memory reallocation
	// This is for Ct
	std::fill( Matrix_elements.begin(), Matrix_elements.begin() + (numFacetsTot+numCells-1), 0);	
    std::fill( Matrix_elements.begin() + numFacetsTot+numCells, Matrix_elements.end(), 2);
    
    for (auto& facet_it : this->M_mesh.getFractureFacetsIdsVector()){
		
		//This is for -T
		auto f_neighbors = facet_it.getFractureNeighbors();
		UInt N_neighbors = 0;
        auto sum_maker = [&N_neighbors](std::pair< Fracture_Juncture, std::vector<UInt> > p){N_neighbors += p.second.size() ;};
        std::for_each(f_neighbors.begin(), f_neighbors.end(), sum_maker);
		Matrix_elements[ numFacetsTot+numCells + facet_it.getFractureId() ] += N_neighbors + 1;
		
		// This is for C and coupling terms on M
		UInt Id_plus   = facet_it.getId();
		UInt Id_minus  = numFacets + facet_it.getFractureId();
		Matrix_elements[ Id_plus ]   += 2;
		Matrix_elements[ Id_minus ]  += 2;
		
	}
    
    S.reserve(Matrix_elements);
    Matrix_elements.clear();
    Matrix_elements.shrink_to_fit();


    // Assemble fractures with FV and impose coupling conditions
    for (auto& facet_it : this->M_mesh.getFractureFacetsIdsVector()){
		
        auto& F_permeability = M_properties.getProperties(facet_it.getZoneCode()).M_permeability;
        const Real Kfn = F_permeability->operator()(0,0);
        const Real F_aperture = M_properties.getProperties(facet_it.getZoneCode()).M_aperture;
		const Real eta = F_aperture / Kfn;
		
		const UInt Id_plus         = facet_it.getId();
		const UInt Id_F            = facet_it.getFractureId();
		const UInt Id_minus        = numFacets + Id_F;
		const Real facetMeasure    = facet_it.getSize();
		
		// Coupling terms on M
		S.coeffRef(Id_plus,Id_plus)    +=   eta*facetMeasure/4. + eta*xsi0*facetMeasure;
		S.coeffRef(Id_plus,Id_minus)   +=   eta*facetMeasure/4. - eta*xsi0*facetMeasure;
		S.coeffRef(Id_minus,Id_minus)  +=   eta*facetMeasure/4. + eta*xsi0*facetMeasure;
		S.coeffRef(Id_minus,Id_plus)   +=   eta*facetMeasure/4. - eta*xsi0*facetMeasure;
		
		// Coupling matrix Ct and C
		S.insert(Id_plus, numFacetsTot+numCells + Id_F)     =  facetMeasure;
		S.insert(Id_minus, numFacetsTot+numCells + Id_F)    = -facetMeasure;
		S.insert(numFacetsTot+numCells + Id_F, Id_plus)     =  facetMeasure;
		S.insert(numFacetsTot+numCells + Id_F, Id_minus)    = -facetMeasure;
	
		// Assemble -T matrix taking into account fractures intersections with the "star-delta" transformation
        for (auto juncture_it : facet_it.getFractureNeighbors()){
			
            std::vector<Real> alphas;
            const Real alphaF = findFracturesAlpha (juncture_it.first, facet_it.getFractureId());

            for (auto neighbors_it : juncture_it.second)
                alphas.emplace_back(findFracturesAlpha (juncture_it.first, neighbors_it));

            Real a_sum = alphaF;
            auto sum_maker = [&a_sum](Real item){a_sum += item;};
            std::for_each(alphas.begin(), alphas.end(), sum_maker);

            for (UInt counter = 0; counter < alphas.size(); ++counter){
				
                const Real QFf = alphaF * alphas[counter] * M_properties.getMobility() / a_sum;
                
                S.coeffRef( numFacetsTot+numCells + Id_F, numFacetsTot+numCells + Id_F ) -= QFf;
                
                S.insert( numFacetsTot+numCells + Id_F, numFacetsTot+numCells + juncture_it.second[counter]) = QFf; 
                
            }
            alphas.clear();
        }
    }
                            
    S.makeCompressed();                      

    // assemble BC on porous matrix
    for (auto& facet_it : this->M_mesh.getBorderFacetsIdsVector()){
		
        const UInt borderId = facet_it.getBorderId();
        const UInt facetId = facet_it.getId();

        if(M_bc.getBordersBCMap().at(borderId).getBCType() == Neumann ){
			
            const Real vel_N = M_bc.getBordersBCMap().at(borderId).getBC()(facet_it.getCentroid());
            rhs[facetId] = vel_N;
        }
        
        else if(M_bc.getBordersBCMap().at(borderId).getBCType() == Dirichlet){
			
			Real pD = M_bc.getBordersBCMap().at(borderId).getBC()(facet_it.getCentroid());
			rhs[facetId] *= pD;
		}
    }

    // assemble BC on fractures
    for (auto& edge_it : this->M_mesh.getBorderTipEdgesIdsVector()){
		
        bool isD = false;
        UInt borderId = 0;

        // select which BC to apply
        for(auto border_it : edge_it.getBorderIds()){
			
            // BC = D > N && the one with greatest id
            if(M_bc.getBordersBCMap().at(border_it).getBCType() == Dirichlet){
				
                if(!isD){
                    isD = true;
                    borderId = border_it;
                }
                else
                    borderId = (border_it > borderId) ? border_it : borderId;
            }
            else if(!isD && M_bc.getBordersBCMap().at(border_it).getBCType() == Neumann){
				
                borderId = (border_it > borderId) ? border_it : borderId;
			}
        }

        // loop over the fracture facets to impose BC on fractures
        for(auto facet_it : edge_it.getSeparatedFacetsIds()){
			 
            if(facetVectorRef[facet_it].isFracture()){
				
                const UInt neighborId = facetVectorRef[facet_it].getFractureFacetId();
                const Real aperture = M_properties.getProperties( facetVectorRef[facet_it].getZoneCode() ).M_aperture;

                if(M_bc.getBordersBCMap().at(borderId).getBCType() == Neumann){
					
                    const Real Q1o = M_bc.getBordersBCMap().at(borderId).getBC()(edge_it.getCentroid()) * edge_it.getSize() * aperture;
                    rhs[ numFacetsTot+numCells + neighborId ] += Q1o;
                } // if
                
                else if(M_bc.getBordersBCMap().at(borderId).getBCType() == Dirichlet){
					
                    const Real alpha1 = findAlpha (facet_it, &edge_it);
                    const Real alpha2 = findDirichletAlpha (facet_it, &edge_it);

                    const Real T12 = 2 * alpha1*alpha2/(alpha1 + alpha2);
                    const Real Q12 = T12 * M_properties.getMobility();
                    const Real Q1o = Q12 * M_bc.getBordersBCMap().at(borderId).getBC()(edge_it.getCentroid());

                    S.coeffRef( numFacetsTot+numCells + neighborId, numFacetsTot+numCells + neighborId) -= Q12; 
                    rhs[ numFacetsTot+numCells + neighborId ] -= Q1o;

                } // else if
            } // if
        }// for
    } // for
    
    std::cout<<" Assembling ended"<<std::endl;
    
    
} // StiffMatrix::assembleMFD

} // namespace FVCode3D
