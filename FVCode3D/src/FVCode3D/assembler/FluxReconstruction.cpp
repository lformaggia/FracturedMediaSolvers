/*!
 * @file FluxReconstruction.cpp
 * @brief This class computes the flux/velocity from the Darcy problem (definitions).
 */

#include <FVCode3D/assembler/FluxReconstruction.hpp>
#include <FVCode3D/property/Properties.hpp>
#include <FVCode3D/property/Permeability.hpp>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/LU>
#include <unsupported/Eigen/SparseExtra>
#include <iostream>

namespace FVCode3D
{

void FluxReconstruction::reconstruct()
{
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
} // FluxReconstruction::reconstruct

} // namespace FVCode3D
