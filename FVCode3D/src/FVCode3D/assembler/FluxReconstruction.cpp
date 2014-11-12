/*!
 * @file FluxReconstruction.cpp
 * @brief This class computes the flux/velocity from the Darcy problem (definitions).
 */

#include <FVCode3D/assembler/FluxReconstruction.hpp>
#include <FVCode3D/property/Properties.hpp>
#include <FVCode3D/property/Permeability.hpp>

#include <FVCode3D/export/ExportVTU.hpp>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/LU>
#include <unsupported/Eigen/SparseExtra>

namespace FVCode3D
{

void FluxReconstruction::localReconstruct(const UInt facetId, const Vector & pressure)
{

}

void FluxReconstruction::reconstruct(const Vector & pressure)
{
    std::cout<<std::endl;
    using Eigen::Dynamic;
    using Eigen::RowMajor;
    using Eigen::ColMajor;

    std::vector<Triplet> Matrix_elements;
    std::vector<UInt> ZMatrix_elements;
    std::vector<UInt> BMatrix_elements;
    Real tCoeff=2.;

    // Local Matrices for Mimetic
    Eigen::Matrix<Real,Dynamic,3> Np;         // facet normals
    Eigen::Matrix<Real,Dynamic,3> Rp;         // centroid * area
    Eigen::Matrix<Real,Dynamic,Dynamic> Z0p;  // Component of Z matrix
    Eigen::Matrix<Real,Dynamic,Dynamic> Z1p;  // Component of Z Matrix
    Eigen::Matrix<Real,Dynamic,Dynamic> Zp;   // Z matrix for internal product= \f$ M^{-1}\f$
    Eigen::Matrix<Real,1,Dynamic> Bpmod;      // areas, all positive
    Eigen::Matrix<Real,Dynamic,Dynamic> Qp;   // Base for the column space of Rp

    // Extract info of mesh. auto&& resolves const &
    auto&& cellVectorRef  = this->M_mesh.getCellsVector();
    auto&& facetVectorRef = this->M_mesh.getFacetsVector();
    const UInt numFacets = facetVectorRef.size();
    const UInt numCells = cellVectorRef.size();
    const UInt numCellsTot = cellVectorRef.size() + this->M_mesh.getFractureFacetsIdsVector().size();

    // Sizing global matrices
    Eigen::SparseMatrix<Real, RowMajor> Z;  // Z matrix for internal product= \f$ M^{-1}\f$
    Eigen::SparseMatrix<Real, RowMajor> Bmod;  // B matrix for signed area
    Eigen::SparseMatrix<Real, RowMajor> Bdmod; // B matrix for Dirichlet area
    Eigen::SparseMatrix<Real, ColMajor> T(numCells,numCells);
    Eigen::SparseMatrix<Real, ColMajor> Td(numCells,numFacets); // T matrix that contains the Dirichlet contributes
    ZMatrix_elements.resize(numFacets,0);
    BMatrix_elements.resize(numCells,0);
    Z.resize(numFacets,numFacets);
    Bmod.resize(numCells,numFacets);
    Bdmod.resize(numFacets,numFacets);

    // First loop to size matrices and avoid memory realloc
    for (auto&& cell : cellVectorRef)
    {
        const std::vector<UInt> & cellFacetsId( cell.getFacetsIds() );
        UInt numCellFacets = cellFacetsId.size();
        UInt cellId        = cell.getId();

        for(UInt localFacetId=0; localFacetId<numCellFacets; ++localFacetId)
        {
            UInt globalFacetId = cellFacetsId[localFacetId];
            const Rigid_Mesh::Facet & fac = facetVectorRef[globalFacetId];
            if(
                ( fac.isBorderFacet() &&
                  (M_bc.getBordersBCMap().at(fac.getBorderId()).getBCType() == Neumann)
                ) ||
                ( fac.isFracture() )
              )
            {
                continue;
            }
            BMatrix_elements[cellId] += 1;
#ifdef DIAGONALZ
            ZMatrix_elements[globalFacetId] += 1;
#else // DIAGONALZ
            ZMatrix_elements[globalFacetId] += numCellFacets;
#endif // DIAGONALZ
        }
    }
    Z.reserve(ZMatrix_elements);
    Bmod.reserve(BMatrix_elements);
    ZMatrix_elements.clear();
    ZMatrix_elements.shrink_to_fit();
    BMatrix_elements.clear();
    BMatrix_elements.shrink_to_fit();

    // Loop on cells
    for (auto&& cell : cellVectorRef)
    {
        const std::vector<UInt> & cellFacetsId( cell.getFacetsIds() );
        UInt numCellFacets  = cellFacetsId.size();
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
        Bpmod.resize(Eigen::NoChange,numCellFacets);
        Bpmod.setZero();

        for(UInt localFacetId=0; localFacetId<numCellFacets; ++localFacetId)
        {
            UInt globalFacetId = cellFacetsId[localFacetId];
            const Rigid_Mesh::Facet & fac = facetVectorRef[globalFacetId];

            if(
                ( fac.isBorderFacet() &&
                  (M_bc.getBordersBCMap().at(fac.getBorderId()).getBCType() == Neumann)
                ) ||
                ( fac.isFracture() )
              )
            {
                continue;
            }

            Point3D facetBaricenter = fac.getCentroid();
            Point3D facetNormal     = fac.getUnsignedNormal();
            Point3D g   = facetBaricenter - cellBaricenter;
            Real facetMeasure       = fac.area();
            Real dotp   = dotProduct(g,facetNormal);
            Real alpha  = (dotp >=0. ? 1.0 : -1.0);

            //! @todo I use the formulation of Nicola (to be reviewed)
            // BEWARE FOR THE CONDITION ON VELOCITY I NEED THE AVERAGE
            // VELOCITY NOT THE FLUX
            Np(localFacetId,0)    = facetNormal[0];
            Np(localFacetId,1)    = facetNormal[1];
            Np(localFacetId,2)    = facetNormal[2];
            Bpmod(0,localFacetId) = facetMeasure;

            g[0] = std::fabs(g[0]) < 1e-15 ? 0 : g[0];
            g[1] = std::fabs(g[1]) < 1e-15 ? 0 : g[1];
            g[2] = std::fabs(g[2]) < 1e-15 ? 0 : g[2];

            Rp(localFacetId,0) = alpha * g[0] * facetMeasure;
            Rp(localFacetId,1) = alpha * g[1] * facetMeasure;
            Rp(localFacetId,2) = alpha * g[2] * facetMeasure;

            if (fac.isBorderFacet() && (M_bc.getBordersBCMap().at(fac.getBorderId()).getBCType() == Dirichlet))
            {
                Bdmod.coeffRef(globalFacetId,globalFacetId) = facetMeasure;
            }
        }

        Qp = Rp.fullPivLu().image(Rp);

        Eigen::Matrix<Real,Dynamic,Dynamic> RtR = Qp.transpose() * Qp;
        Eigen::Matrix<Real,Dynamic,Dynamic> RRtRiRt = Qp * RtR.inverse() * Qp.transpose();

/*
        std::stringstream  ss;
        ss << "RRtRiRt_" << cellId << ".m";
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
        ss << "QQt_" << cellId << ".m";
        Eigen::saveMarket( QQt, ss.str() );
        ss.str(std::string());
        ss << "Rp_" << cellId << ".m";
        Eigen::saveMarket( Rp, ss.str() );
        ss.str(std::string());
        ss << "Qp_" << cellId << ".m";
        Eigen::saveMarket( Qp, ss.str() );
*/

        Z0p = ( K->norm() * ( 1. / cellMeasure ) ) *
              (Np * Np.transpose());
        Z1p = ( tCoeff * ( K->operator()(0,0) + K->operator()(1,1) + K->operator()(2,2) ) *
              ( 1. / cellMeasure ) ) *
              ( Eigen::MatrixXd::Identity(numCellFacets,numCellFacets) -
                RRtRiRt
//                Qp * Qp.transpose()
              );

        Zp  = Z0p + Z1p;

        for(UInt iloc=0; iloc<numCellFacets; ++iloc)
        {
            const UInt i = cellFacetsId[iloc];
            const Real Bcoeffmod = Bpmod(0,iloc);

            if(Bcoeffmod != 0.)
            {
                Bmod.insert(cellId,i) = Bcoeffmod;
            }

            Real Zcoeff = Zp(iloc,iloc);

            if( Zcoeff != 0. )
            {
                Z.coeffRef(i,i) += Zp(iloc,iloc);
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
        }
    }

    Real Df;
    Real T1f, Q1f, T2f, Q2f, QFf, Q1o;
    Real alpha1, alpha2, alphaF, alphaf;
    UInt neighbor1id, neighbor2id;
    std::vector<Real> alphas;

    // assemble fractures with FV
    for (auto facet_it : this->M_mesh.getFractureFacetsIdsVector())
    {
        auto& F_permeability = M_properties.getProperties(facet_it.getZoneCode()).M_permeability;
        Real F_aperture = M_properties.getProperties(facet_it.getZoneCode()).M_aperture;
        Df = F_aperture/2.;
        alphaf = facet_it.getSize() * F_permeability->norm() * ( 1. / Df );

        neighbor1id = facet_it.getSeparatedCellsIds()[0];
        neighbor2id = facet_it.getSeparatedCellsIds()[1];

        alpha1 = findAlpha(neighbor1id, &facet_it);
        alpha2 = findAlpha(neighbor2id, &facet_it);

        T1f = alpha1*alphaf/(alpha1 + alphaf);
        Q1f = T1f * M_properties.getMobility();

        T2f = alphaf*alpha2/(alphaf + alpha2);
        Q2f = T2f * M_properties.getMobility();

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
            alphaF = findFracturesAlpha (juncture_it.first, facet_it.getFractureId());

            for (auto neighbors_it : juncture_it.second)
            {
                alphas.emplace_back(findFracturesAlpha (juncture_it.first, neighbors_it));
            }

            Real a_sum = alphaF;
            auto sum_maker = [&a_sum](Real item){a_sum += item;};
            std::for_each(alphas.begin(), alphas.end(), sum_maker);

            for (UInt counter = 0; counter < alphas.size(); ++counter)
            {
                QFf = alphaF*alphas[counter] * M_properties.getMobility() / a_sum;
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
    for (auto facet_it : this->M_mesh.getBorderFacetsIdsVector())
    {
        neighbor1id = facet_it.getSeparatedCellsIds()[0];
        UInt borderId = facet_it.getBorderId();
        UInt facetId = facet_it.getId();

        if(M_bc.getBordersBCMap().at(borderId).getBCType() == Neumann )
        {
            // Nella cond di bordo c'è già il contributo della permeabilità e mobilità (ma non la densità!)
            Q1o = M_bc.getBordersBCMap().at(borderId).getBC()(facet_it.getCentroid()) * facet_it.getSize();

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

                    const Real T12 = alpha1*alpha2/(alpha1 + alpha2);
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
    Eigen::saveMarket( *(this->M_Matrix), "A.m" );
    Eigen::saveMarket( *M_b, "RHS.m" );
    //this->M_Matrix->setFromTriplets(Matrix_elements.begin(), Matrix_elements.end());
    std::cout<<" Assembling ended"<<std::endl;










    // Extract info of mesh. auto&& resolves const &
    auto&& fractureVectorRef = this->M_mesh.getFractureFacetsIdsVector();
    auto&& borderVectorRef = this->M_mesh.getBorderFacetsIdsVector();
    const UInt numFractures = fractureVectorRef.size();
    const UInt numFacetsTot = numFacets + numFractures;

    // Resize the vectors
    M_flux = Vector::Constant(numFacetsTot, 0.);
    //M_velocity.resize(3*numFacetsTot);

    M_flux = - Z * Bmod.transpose() * pressure;

    Vector pressureD(Vector::Constant(Bdmod.cols(), 0.));
    for (auto facet_it : borderVectorRef)
    {
        UInt borderId = facet_it.getBorderId();
        UInt facetId = facet_it.getId();

        if(M_bc.getBordersBCMap().at(borderId).getBCType() == Neumann )
        {
            // Nella cond di bordo c'è già il contributo della permeabilità e mobilità (ma non la densità!)
            M_flux(facetId) = M_bc.getBordersBCMap().at(borderId).getBC()(facet_it.getCentroid()) * facet_it.getSize();
        }
        else if(M_bc.getBordersBCMap().at(borderId).getBCType() == Dirichlet)
        {
            pressureD(facetId) = M_bc.getBordersBCMap().at(borderId).getBC()(facet_it.getCentroid());
        }
    }

    M_flux += - Z * Bdmod.transpose() * pressureD;

    for(UInt i=0; i<(numFacets + numFractures); ++i)
    {
        M_flux(i) /= 2;
    }

    ExporterVTU exporter;
    exporter.exportFlux(this->M_mesh, "./flux.vtu", M_flux, "flux");
} // FluxReconstruction::reconstruct


//    using Eigen::Dynamic;
//    using Eigen::RowMajor;
//    using Eigen::ColMajor;
//
//    // Extract info of mesh. auto&& resolves const &
//    auto&& cellVectorRef   = this->M_mesh.getCellsVector();
//    auto&& facetVectorRef  = this->M_mesh.getFacetsVector();
//    UInt numFacets = facetVectorRef.size();
//    UInt numCells = cellVectorRef.size();
//    UInt numFractures = this->M_mesh.getFractureFacetsIdsVector().size();
//    UInt numCellsTot = cellVectorRef.size() + numFractures;
//
//    // Resize the vectors
//    M_flux.resize(numFacets + numFractures);
//    M_velocity.resize(3*numCellsTot);
//
//    Real tCoeff=2.;
//
//    std::vector<Triplet> Matrix_elements;
//    Real alpha1, alpha2, alphaF, alphaf;
//    UInt neighbor1id, neighbor2id;
//    std::vector<Real> alphas;
//    Point3D f1, f2;
//    Real Df;
//    Real T12, Q12, T1f, Q1f, T2f, Q2f, QFf, Q1o;
//
//    // Local Matrices for Mimetic
//    Eigen::Matrix<double,Dynamic,3> Np; // facet normals
//    Eigen::Matrix<double,Dynamic,Dynamic> Z0p;// COmpunent of Z matrix
//    Eigen::Matrix<double,Dynamic,Dynamic> Z1p;// Component of Z Matrix
//    Eigen::Matrix<double,Dynamic,Dynamic> Zp; //! Z matrix for internal product= \f$ M^{-1}\f$
//    Eigen::DiagonalMatrix<double,Dynamic> Bp; // B Matrix
//    Eigen::Matrix<double,Dynamic,3> Rp; // Matrix R
//    Eigen::Matrix<double,Dynamic,Dynamic> Qp; // base for the column space of Rp
//    Eigen::Matrix<double,Dynamic,Dynamic> Tp; // T matrix
//    Eigen::Matrix<double,Dynamic,1> Pp; // pressure local matrix
//    Eigen::Matrix<double,Dynamic,1> Fp; // fluxes local matrix
//
//    std::map<UInt,UInt> facetToCell;
//
//    // Loop on cells
//    std::cout<<" Starting loops on cells"<<std::endl;
//    for (auto&& cell : cellVectorRef)
//    {
//        std::vector<UInt> const & cellFacetsId( cell.getFacetsIds() );
//        UInt numCellFacets = cellFacetsId.size();
//        UInt numLocalCells = numCellFacets+1;
//        auto K = M_properties.getProperties(cell.getZoneCode()).M_permeability;
//        auto cellBaricenter = cell.getCentroid();
//        auto cellMeasure    = cell.getVolume();
//        auto cellId         = cell.getId();
//        if(cellId % 500 == 0)
//        {
//            std::cout<<"Done "<< cellId<<" Cells"<<std::endl;
//        }
//        // Resize local matrices
//        Np.resize(numCellFacets,Eigen::NoChange);
//        Np.setZero();
//        Rp.resize(numCellFacets,Eigen::NoChange);
//        Rp.setZero();
//        Zp.resize(numCellFacets,numCellFacets);
//        Zp.setZero();
//        Z0p.resize(numCellFacets,numCellFacets);
//        Z0p.setZero();
//        Z1p.resize(numCellFacets,numCellFacets);
//        Z1p.setZero();
//        Bp.resize(numCellFacets);
//        Bp.setZero();
//        Pp.resize(numCellFacets,Eigen::NoChange);
//        Pp.setZero();
//        Fp.resize(numCellFacets,Eigen::NoChange);
//        Fp.setZero();
//
//        for(UInt localFacetId=0; localFacetId<numCellFacets; ++localFacetId)
//        {
//            UInt globalFacetId = cellFacetsId[localFacetId];
//            const Rigid_Mesh::Facet & fac=facetVectorRef[globalFacetId];
//            if(fac.isBorderFacet() && (M_bc.getBordersBCMap().at(fac.getBorderId()).getBCType() == Neumann))
//            {
//                continue;
//            }
//            Real alpha(0.);
//            Point3D facetBaricenter= fac.getCentroid();
//            Point3D facetNormal    = fac.getUnsignedNormal();
//            auto facetMeasure   = fac.area();
//            Point3D g = facetBaricenter-cellBaricenter;
//            Real dotp = dotProduct(g,facetNormal);
//            alpha = (dotp >=0.? 1.0:-1.0);
//
//            UInt oppositeCell;
//            if(fac.isBorderFacet())
//            {
//                oppositeCell = numCellsTot + fac.getBorderId();
//            }
//            else
//            {
//                oppositeCell = fac.getSeparatedCellsIds()[0] != cellId ? fac.getSeparatedCellsIds()[0] :
//                fac.getSeparatedCellsIds()[1];
//            }
//            facetToCell.insert(std::pair<UInt,UInt>(localFacetId, oppositeCell));
//
//            Np(localFacetId,0)=facetNormal[0];
//            Np(localFacetId,1)=facetNormal[1];
//            Np(localFacetId,2)=facetNormal[2];
//            Bp.diagonal()[localFacetId]=alpha*facetMeasure;
//
//            g[0] = std::fabs(g[0]) < 1e-15 ? 0 : g[0];
//            g[1] = std::fabs(g[1]) < 1e-15 ? 0 : g[1];
//            g[2] = std::fabs(g[2]) < 1e-15 ? 0 : g[2];
//
//            Rp(localFacetId,0)=alpha*g[0]*facetMeasure;
//            Rp(localFacetId,1)=alpha*g[1]*facetMeasure;
//            Rp(localFacetId,2)=alpha*g[2]*facetMeasure;
//            if (fac.isBorderFacet() && (M_bc.getBordersBCMap().at(fac.getBorderId()).getBCType() == Dirichlet))
//            {
//                Bd.coeffRef(globalFacetId,globalFacetId) = alpha*facetMeasure;
//                count++;
//            }
//        }
//
//        Qp = Rp.fullPivLu().image(Rp);
//
//        Z0p = (K * ( 1. / cellMeasure ) ) *
//              (Np * Np.transpose());
//        Z1p = (tCoeff * K * ( 1. / cellMeasure) ) *
//                ( Eigen::MatrixXd::Identity(numCellFacets,numCellFacets) -
//                Qp * Qp.transpose() );
//
//        Zp  = Z0p + Z1p;
//
//        Tp = - Bp * Zp * Bp.transpose();
//        for(UInt localFacetId=0; localFacetId<numCellFacets;++localFacetId)
//        {
//            if(facetToCell[localFacetId] < numCellsTot)
//            {
//                Pp(localFacetId,0) = M_pressure[facetToCell[localFacetId]];
//            }
//            else
//            {
//                Pp(localFacetId,0) = M_bc.getBordersBCMap().at(facetToCell[localFacetId]-numCellsTot).getBC()(facetVectorRef[cellFacetsId[localFacetId]].getCentroid());
//            }
//        }
//        Fp = Tp * Pp;
//
//        for(UInt localFacetId=0; localFacetId<numCellFacets; ++localFacetId)
//        {
//            M_flux[facetToCell[localFacetId]] = Pp(localFacetId) > Fp(localFacetId)
//        }
//    }
//

} // namespace FVCode3D
