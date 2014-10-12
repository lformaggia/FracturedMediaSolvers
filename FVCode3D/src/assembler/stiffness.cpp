/*!
 *  @file stiffness.cpp
 *  @brief This class build a Stiffness-matrix of the Darcy problem (definitions).
 */

#include "assembler/stiffness.hpp"
#include "property/Properties.hpp"
#include "Eigen/Core"
#include "Eigen/SparseCore"
#include "geometry/Point3D.hpp"
#include <Eigen/LU>
#include <iostream>

#include <unsupported/Eigen/SparseExtra>

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
    using Eigen::Dynamic;
    using Eigen::RowMajor;
    using Eigen::ColMajor;
    using Geometry::Rigid_Mesh;

    std::vector<Triplet> Matrix_elements;
    std::vector<UInt> ZMatrix_elements;
    std::vector<UInt> BMatrix_elements;
    Real tCoeff=2.;
    Real alpha1, alpha2, alphaF, alphaf;
    UInt neighbor1id, neighbor2id;
    std::vector<Real> alphas;
    Generic_Vector f1, f2;
    Real Df;
    Real T12, Q12, T1f, Q1f, T2f, Q2f, QFf, Q1o;

    // Local Matrices for Mimetic
    Eigen::Matrix<double,Dynamic,3> Np; // facet normals
    Eigen::Matrix<double,Dynamic,Dynamic> Z0p;// COmpunent of Z matrix
    Eigen::Matrix<double,Dynamic,Dynamic> Z1p;// Component of Z Matrix
    Eigen::Matrix<double,Dynamic,Dynamic> Zp; //! Z matrix for internal product= \f$ M^{-1}\f$
    Eigen::Matrix<double,1,Dynamic> Bp; // B Matrix
    Eigen::Matrix<double,1,Dynamic> Bpmod; // B Matrix, all positive
    Eigen::Matrix<double,Dynamic,3> Rp; // Matrix R
    Eigen::Matrix<double,Dynamic,Dynamic> Q; // Base for the column space of Rp

    // Extract info of mesh. auto&& resolves const &
    auto&& cellVectorRef   = this->M_mesh.getCellsVector();
    auto&& facetVectorRef  = this->M_mesh.getFacetsVector();
    UInt numFacets = facetVectorRef.size();
    UInt numCells = cellVectorRef.size();
    UInt numCellsTot = cellVectorRef.size() + this->M_mesh.getFractureFacetsIdsVector().size();

    // SIZING GLOBAL MATRICES
    //Eigen::SparseMatrix<double,RowMajor> Z(numFacets,numFacets); //! Z matrix for internal product= \f$ M^{-1}\f$
    //Eigen::SparseMatrix<double,RowMajor> Bd(numFacets,numFacets);// Dirichlet B matrix
    //Eigen::SparseMatrix<double,RowMajor> B(numCells,numFacets);
    Eigen::SparseMatrix<double,ColMajor> Td(numCells,numFacets); // T matrix that contains the Dirichlet contributes
    Eigen::SparseMatrix<double,ColMajor> T(numCells,numCells);
    ZMatrix_elements.resize(numFacets,0);
    BMatrix_elements.resize(numCells,0);
    Z.resize(numFacets,numFacets);
    B.resize(numCells,numFacets);
    Bmod.resize(numCells,numFacets);
    Bd.resize(numFacets,numFacets);
    Bdmod.resize(numFacets,numFacets);

    // First loop to size matrices and avoid memory realloc
    for (auto&& cell : cellVectorRef)
    {
        std::vector<UInt> const & cellFacetsId( cell.getFacetsIds() );
        UInt numCellFacets = cellFacetsId.size();
        UInt cellId         = cell.getId();

        for(UInt localFacetId=0; localFacetId<numCellFacets;++localFacetId)
        {
            UInt globalFacetId = cellFacetsId[localFacetId];
            Rigid_Mesh::Facet const & fac=facetVectorRef[globalFacetId];
            if(fac.isBorderFacet() && m_Bc.getBordersBCMap().at(fac.getBorderId()).getBCType() == Neumann)
                continue;
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
    Bmod.reserve(BMatrix_elements);
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
        auto K = M_properties.getProperties(cell.getZoneCode()).M_permeability;
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
        Bpmod.resize(Eigen::NoChange,numCellFacets);
        Bpmod.setZero();
        for(UInt localFacetId=0; localFacetId<numCellFacets;++localFacetId)
        {

            UInt globalFacetId = cellFacetsId[localFacetId];
            Rigid_Mesh::Facet const & fac=facetVectorRef[globalFacetId];
            if(fac.isBorderFacet() && (m_Bc.getBordersBCMap().at(fac.getBorderId()).getBCType() == Neumann))
                continue;
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
            Np(localFacetId,2)=facetNormal[2];
            Bp(0,localFacetId)=alpha*facetMeasure;
            Bpmod(0,localFacetId)=facetMeasure;

            g[0] = std::fabs(g[0]) < 1e-15 ? 0 : g[0];
            g[1] = std::fabs(g[1]) < 1e-15 ? 0 : g[1];
            g[2] = std::fabs(g[2]) < 1e-15 ? 0 : g[2];

            Rp(localFacetId,0)=alpha*g[0]*facetMeasure;
            Rp(localFacetId,1)=alpha*g[1]*facetMeasure;
            Rp(localFacetId,2)=alpha*g[2]*facetMeasure;
            if (fac.isBorderFacet() && (m_Bc.getBordersBCMap().at(fac.getBorderId()).getBCType() == Dirichlet))
            {
                Bd.coeffRef(globalFacetId,globalFacetId) = alpha*facetMeasure;
                Bdmod.coeffRef(globalFacetId,globalFacetId) = facetMeasure;
            }
        }

        Q = Rp.fullPivLu().image(Rp);

        Z0p = (K/cellMeasure)*(Np*Np.transpose());
        Z1p = (tCoeff*K/cellMeasure)*(
                Eigen::MatrixXd::Identity(numCellFacets,numCellFacets)-
                Q*Q.transpose());

        Zp  = Z0p + Z1p;

        for(UInt iloc=0; iloc<numCellFacets;++iloc)
        {
            UInt i=cellFacetsId[iloc];
            Real Bcoeff=Bp(0,iloc);
            Real Bcoeffmod=Bpmod(0,iloc);
            if(Bcoeff != 0.)
            {
                B.insert(cellId,i)=Bcoeff;
                Bmod.insert(cellId,i)=Bcoeffmod;
            }
            Real Zcoeff = Zp(iloc,iloc);
            if(Zcoeff != 0. )
                Z.coeffRef(i,i)+=Zp(iloc,iloc);
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
    std::cout<<"Building Trasmissibility Dirichelt matrix (global)"<<std::endl;
    Td = B * Z * Bd.transpose();
    Td.makeCompressed();
    // Clean up to seve memory (I hope)
    std::cout<<" Matrix T has "<<T.rows()<<" rows and "<<T.cols()<<" columns and "<<T.nonZeros()<<" Non zeros"<<std::endl;
    std::cout<<" Num Cells"<<numCells<<std::endl;
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
        UInt facetId = facet_it.getFacetId();

        if(m_Bc.getBordersBCMap().at(borderId).getBCType() == Neumann )
        {
            // Nella cond di bordo c'è già il contributo della permeabilità e mobilità (ma non la densità!)
            Q1o = m_Bc.getBordersBCMap().at(borderId).getBC()(facet_it.getCentroid()) * facet_it.getSize();

            _b->operator()(neighbor1id) += Q1o;
        }
        else if(m_Bc.getBordersBCMap().at(borderId).getBCType() == Dirichlet)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(Td,facetId); it; ++it)
                _b->operator()(it.row()) += it.value() * m_Bc.getBordersBCMap().at(borderId).getBC()(facet_it.getCentroid());
        }
    }

    // Put everything together
    Eigen::SparseMatrix<double,ColMajor> Mborder(numCellsTot,numCellsTot);
    Mborder.setFromTriplets(Matrix_elements.begin(), Matrix_elements.end());
    this->M_Matrix->resize(numCellsTot,numCellsTot);
    T.conservativeResize(numCellsTot,numCellsTot);
    *(this->M_Matrix) = T + Mborder;
    //this->M_Matrix->setFromTriplets(Matrix_elements.begin(), Matrix_elements.end());
    std::cout<<" Assembling ended"<<std::endl;
}

void StiffMatrix::reconstructFlux(const Vector & pressure)
{
    using Eigen::Dynamic;
    using Eigen::RowMajor;
    using Eigen::ColMajor;
    using Geometry::Rigid_Mesh;

    // Extract info of mesh. auto&& resolves const &
    auto&& facetVectorRef  = this->M_mesh.getFacetsVector();
    auto&& fractureVectorRef = this->M_mesh.getFractureFacetsIdsVector();
    auto&& borderVectorRef = this->M_mesh.getBorderFacetsIdsVector();
    UInt numFacets = facetVectorRef.size();
    UInt numFractures = fractureVectorRef.size();
    UInt numFacetsTot = numFacets + numFractures;

    // Resize the vectors
    M_flux = Vector::Constant(numFacetsTot, 0.);
    //M_velocity.resize(3*numFacetsTot);

    M_flux = - Z * Bmod.transpose() * pressure;

    Vector pressureD(Vector::Constant(Bdmod.cols(), 0.));
    for (auto facet_it : borderVectorRef)
    {
        UInt neighbor1id = facet_it.getSeparated()[0];
        UInt borderId = facet_it.getBorderId();
        UInt facetId = facet_it.getFacetId();

        if(m_Bc.getBordersBCMap().at(borderId).getBCType() == Neumann )
        {
            // Nella cond di bordo c'è già il contributo della permeabilità e mobilità (ma non la densità!)
            M_flux(facetId) = m_Bc.getBordersBCMap().at(borderId).getBC()(facet_it.getCentroid()) * facet_it.getSize();
        }
        else if(m_Bc.getBordersBCMap().at(borderId).getBCType() == Dirichlet)
        {
            pressureD(facetId) = m_Bc.getBordersBCMap().at(borderId).getBC()(facet_it.getCentroid());
        }
    }

    M_flux += - Z * Bdmod.transpose() * pressureD;

//     for(UInt i=0; i<(numFacets + numFractures); ++i)
//         M_flux(i) /= 2;

}

} // namespace Darcy
