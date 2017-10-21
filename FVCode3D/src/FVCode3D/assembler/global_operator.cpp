/*!
 * @file global_operator.cpp
 * @brief These classes implement global operators (definitions).
 */

#include <vector>
#include <cmath>
#include <FVCode3D/assembler/global_operator.hpp>
#include <FVCode3D/assembler/local_operator.hpp>
#include <FVCode3D/boundaryCondition/BC.hpp>
#include <FVCode3D/mesh/RigidMesh.hpp>
#include <FVCode3D/property/Permeability.hpp>
#include <Eigen/LU>
#include <unsupported/Eigen/SparseExtra>

namespace FVCode3D
{
	
void global_InnerProduct::reserve_space()      
{
	auto & M = *M_matrix;
	std::vector<UInt> Matrix_elements(Ncol,0);
	for (auto& cell : M_mesh.getCellsVector())
	{
	const std::vector<UInt> & cellFacetsId( cell.getFacetsIds() );
	UInt numCellFacets = cellFacetsId.size();

	for(UInt localFacetId=0; localFacetId<numCellFacets; ++localFacetId)
    {		
		UInt globalFacetId = cellFacetsId[localFacetId];
		const Rigid_Mesh::Facet & fac = M_mesh.getFacetsVector()[globalFacetId];
			
		// This is to take into account the decoupling of fractures facets
		if( fac.isFracture() && (cell.orientationFacet(fac)<0) )
			globalFacetId = M_mesh.getFacetsVector().size() + fac.getFractureFacetId();

		Matrix_elements[globalFacetId] += numCellFacets;        
        }
    }
    M.reserve(Matrix_elements);		
}

void global_InnerProduct::assemble()   
{
	auto & M                = *M_matrix;
	auto & facetVectorRef   = M_mesh.getFacetsVector();
	const UInt numFacets    = facetVectorRef.size();
	
	std::cout<<"Assembling mimetic inner product..."<<std::endl;
	
	for (auto& cell : M_mesh.getCellsVector())
	{	 	
		// Definisco prodotto interno locale
		local_InnerProduct   localIP(M_mesh,cell);
		// Assemblo prodotto interno locale
		localIP.assemble();
		
		auto & Mp = localIP.getMp();
		
		const std::vector<UInt> & cellFacetsId( cell.getFacetsIds() );
		const UInt numCellFacets  = cellFacetsId.size();
		
		// Assemblo matrice locale in matrice globale
		for(UInt iloc=0; iloc<numCellFacets; ++iloc){
			
            UInt i = cellFacetsId[iloc];                                 //global Id
            const Rigid_Mesh::Facet & fac = facetVectorRef[i];     
            
            // This is to take into account the decoupling of fractures facets
			if( fac.isFracture() && (cell.orientationFacet(fac)<0) )
					i = numFacets + fac.getFractureFacetId();

            Real Mcoeff = Mp(iloc,iloc);

            if( Mcoeff != 0. )
                M.coeffRef(i,i) += Mcoeff;

            for(UInt jloc=iloc+1; jloc<numCellFacets; ++jloc){

                UInt j = cellFacetsId[jloc];
                const Rigid_Mesh::Facet & fac = facetVectorRef[j];
                Mcoeff = Mp(iloc,jloc);
                
                // This is to take into account the decoupling of fractures facets
				if( fac.isFracture() && (cell.orientationFacet(fac)<0) )
					j = numFacets + fac.getFractureFacetId();

                if (Mcoeff != 0.){
                    M.coeffRef(i,j) += Mcoeff;
                    M.coeffRef(j,i) += Mcoeff;
                } 
            }
        }
	}
	std::cout<<"Done."<<std::endl;
}

void global_InnerProduct::ImposeBC(SpMat & S, Vector & rhs)      
{
	std::cout<<std::endl;
	std::cout << "Zero the Neumann row of S and impose BC on the rhs" << std::endl<<std::endl;
	
	for( auto facet_it : M_mesh.getBorderFacetsIdsVector() )
	{	
		UInt borderId = facet_it.getBorderId();
		UInt facetId  = facet_it.getId();
			
		if( M_bc.getBordersBCMap().at(borderId).getBCType() == Neumann )
		{
			// Impose Neumann BC on the matrix of the system
			S.prune( [facetId]( UInt i, UInt, Real ){ return ( i!= facetId ); });
			S.coeffRef(facetId,facetId) = 1;
			// Impose Neumann BC on the rhs
			const Real vel_N = M_bc.getBordersBCMap().at(borderId).getBC()(facet_it.getCentroid());
            rhs[facetId] = vel_N;
		}
			
		else if(M_bc.getBordersBCMap().at(borderId).getBCType() == Dirichlet)
		{
			// Impose Dirichlet BC on the rhs
			Real pD = M_bc.getBordersBCMap().at(borderId).getBC()(facet_it.getCentroid());
			const Real facetMeasure = facet_it.getSize();
			const Rigid_Mesh::Cell & cell = M_mesh.getCellsVector()[ facet_it.getSeparatedCellsIds()[0] ];
			const Real alpha = cell.orientationFacet( M_mesh.getFacetsVector()[ facetId ] ); 
			rhs[facetId] = - alpha * facetMeasure * pD;
		}
	}
}

void global_InnerProduct::ImposeBC(Vector & rhs)      
{
	std::cout << "Zero the Neumann row of M and impose BC on the rhs" << std::endl<<std::endl;;
	auto & M = *M_matrix;
	
	for( auto facet_it : M_mesh.getBorderFacetsIdsVector() )
	{	
		UInt borderId = facet_it.getBorderId();
		UInt facetId  = facet_it.getId();
			
		if( M_bc.getBordersBCMap().at(borderId).getBCType() == Neumann )
		{
			// Impose Neumann BC on the matrix of the system
			M.prune( [facetId]( UInt i, UInt, Real ){ return ( i!= facetId ); });
			M.coeffRef(facetId,facetId) = 1;
			// Impose Neumann BC on the rhs
			const Real vel_N = M_bc.getBordersBCMap().at(borderId).getBC()(facet_it.getCentroid());
            rhs[facetId] = vel_N;
		}
			
		else if(M_bc.getBordersBCMap().at(borderId).getBCType() == Dirichlet)
		{
			// Impose Dirichlet BC on the rhs
			Real pD = M_bc.getBordersBCMap().at(borderId).getBC()(facet_it.getCentroid());
			const Real facetMeasure = facet_it.getSize();
			const Rigid_Mesh::Cell & cell = M_mesh.getCellsVector()[ facet_it.getSeparatedCellsIds()[0] ];
			const Real alpha = cell.orientationFacet( M_mesh.getFacetsVector()[ facetId ] ); 
			rhs[facetId] = - alpha * facetMeasure * pD;
		}
	}
}


void global_Div::reserve_space()      
{
	auto & B   = *M_matrix;
	auto & Dt  = *Dt_matrix;
	std::vector<UInt> Matrix_elements(Ncol,0);
	std::vector<UInt> DtMatrix_elements(Nrow,0);
	
	for (auto& cell : M_mesh.getCellsVector())
	{
	const std::vector<UInt> & cellFacetsId( cell.getFacetsIds() );
	const UInt cellId = cell.getId();

	for(UInt localFacetId=0; localFacetId<cellFacetsId.size(); ++localFacetId)
    {		
		UInt globalFacetId = cellFacetsId[localFacetId];
		const Rigid_Mesh::Facet & fac = M_mesh.getFacetsVector()[globalFacetId];
			
		// This is to take into account the decoupling of fractures facets
		if( fac.isFracture() && (cell.orientationFacet(fac)<0) )
			globalFacetId = M_mesh.getFacetsVector().size() + fac.getFractureFacetId();

		Matrix_elements[globalFacetId] += 1;           // space for B element
		
		if( fac.isBorderFacet() &&
		( M_bc.getBordersBCMap().at(fac.getBorderId()).getBCType() == Neumann ) )
			continue;
		DtMatrix_elements[cellId] += 1;                // space for Dt element
		
        }
    }
    B.reserve(Matrix_elements);	
    Dt.reserve(DtMatrix_elements);	
}

void global_Div::assemble()   
{
	auto & B                = *M_matrix;
	auto & Dt               = *Dt_matrix;
	auto & facetVectorRef   = M_mesh.getFacetsVector();
	const UInt numFacets    = facetVectorRef.size();
	
	for (auto& cell : M_mesh.getCellsVector())
	{	 	
		// Definisco divergenza locale
		local_Div   localDIV(M_mesh,cell);
		// Assemblo divergenza locale
		localDIV.assemble();
		auto & Bp = localDIV.getBp();
		
		const UInt cellId = cell.getId();
		const std::vector<UInt> & cellFacetsId( cell.getFacetsIds() );
		const UInt numCellFacets  = cellFacetsId.size();
		
		// Assemblo matrice locale in matrice globale
		for(UInt iloc=0; iloc<numCellFacets; ++iloc){
			
            UInt i = cellFacetsId[iloc];                                 //global Id
            const Rigid_Mesh::Facet & fac = facetVectorRef[i]; 
            const Real Bcoeff = Bp[iloc];    
            
            // This is to take into account the decoupling of fractures facets
			if( fac.isFracture() && (cell.orientationFacet(fac)<0) )
					i = numFacets + fac.getFractureFacetId();

            if( Bcoeff != 0. ){
				B.insert(cellId, i) = Bcoeff;
				
				if( fac.isBorderFacet() &&
				( M_bc.getBordersBCMap().at(fac.getBorderId()).getBCType() == Neumann ) )
					continue;
				Dt.insert(i, cellId) = Bcoeff;
			}
        }
	}
}

void CouplingConditions::reserve_space()
{
	auto & C = *M_matrix;
	std::vector<UInt> CMatrix_elements( M.cols(), 0 );
	std::vector<UInt> MMatrix_elements( M.cols(), 0 );
	
	for (auto& facet_it : M_mesh.getFractureFacetsIdsVector()){
		UInt Id_plus   = facet_it.getId();
		UInt Id_minus  = M_mesh.getFacetsVector().size() + facet_it.getFractureId();
		// This is for C
		CMatrix_elements[ Id_plus ]   = 1;
		CMatrix_elements[ Id_minus ]  = 1;
		// This is for M
		MMatrix_elements[ Id_plus ]   = 1;
		MMatrix_elements[ Id_minus ]  = 1;
	}
	C.reserve(CMatrix_elements);
	M.reserve(MMatrix_elements);
}

void CouplingConditions::assemble()
{
	auto & C = *M_matrix;
    for (auto& facet_it : M_mesh.getFractureFacetsIdsVector())
	{		
        auto& F_permeability = M_mesh.getPropertiesMap().getProperties(facet_it.getZoneCode()).M_permeability;
        const Real Kfn = F_permeability->operator()(0,0);
        const Real F_aperture = M_mesh.getPropertiesMap().getProperties(facet_it.getZoneCode()).M_aperture;
        
		const Real eta          = F_aperture / Kfn;
		const Real xsi0     =  (2. * xsi - 1.) / 4.;
		
		const UInt Id_plus         = facet_it.getId();
		const UInt Id_F            = facet_it.getFractureId();
		const UInt Id_minus        = M_mesh.getFacetsVector().size() + Id_F;
		const Real facetMeasure    = facet_it.getSize();
		
		// Coupling terms on M
		M.coeffRef(Id_plus,Id_plus)    +=   eta*facetMeasure/4. + eta*xsi0*facetMeasure;
		M.insert(Id_plus,Id_minus)      =   eta*facetMeasure/4. - eta*xsi0*facetMeasure;
		M.coeffRef(Id_minus,Id_minus)  +=   eta*facetMeasure/4. + eta*xsi0*facetMeasure;
		M.insert(Id_minus,Id_plus)      =   eta*facetMeasure/4. - eta*xsi0*facetMeasure;
		
		// Coupling matrix C
		C.insert(Id_F, Id_plus)     =  facetMeasure;
		C.insert(Id_F, Id_minus)    = -facetMeasure;
	}
}


Real FluxOperator::findFracturesAlpha (const Fracture_Juncture fj, const UInt n_Id) const
{
    const Point3D borderCenter = getBorderCenter(fj);
    const Point3D cellCenter = M_mesh.getFractureFacetsIdsVector()[n_Id].getCentroid();

    const Real A = M_mesh.getPropertiesMap().getProperties(M_mesh.getFractureFacetsIdsVector()[n_Id].getZoneCode()).M_aperture *
                  (M_mesh.getNodesVector()[fj.second]-M_mesh.getNodesVector()[fj.first]).norm();

    const PermPtr_Type & k =
            M_mesh.getPropertiesMap().getProperties(M_mesh.getFractureFacetsIdsVector()[n_Id].getZoneCode()).M_permeability;

    Point3D f = borderCenter - cellCenter;
    const Real D = sqrt(dotProduct(f, f));
    assert( D != 0 );

    f /= D;

    Point3D normal = crossProduct(f, M_mesh.getNodesVector()[fj.second]-M_mesh.getNodesVector()[fj.first]); // k = f x l
    normal = crossProduct(M_mesh.getNodesVector()[fj.second]-M_mesh.getNodesVector()[fj.first], normal); // n = l x k
    normal.normalize();

    const Real KnDotF = fabs(dotProduct(k * normal, f)); // K n * f

    return A * KnDotF / D;
}

Real FluxOperator::findAlpha (const UInt & facetId, const Edge_ID * edge) const
{
    const Point3D edgeCenter = edge->getCentroid();
    const Point3D facetCenter = M_mesh.getFacetsVector()[facetId].getCentroid();

    const Real A = edge->getSize() *
                   M_mesh.getPropertiesMap().getProperties(M_mesh.getFacetsVector()[facetId].getZoneCode()).M_aperture;

    const PermPtr_Type & k =
            M_mesh.getPropertiesMap().getProperties(M_mesh.getFacetsVector()[facetId].getZoneCode()).M_permeability;

    Point3D f = edgeCenter - facetCenter;
    const Real D = sqrt(dotProduct(f, f));
    assert( D != 0 );

    f /= D;

    Point3D normal = crossProduct(  f,
                            M_mesh.getNodesVector()[edge->getEdge().getVerticesIds()[0]] -
                            M_mesh.getNodesVector()[edge->getEdge().getVerticesIds()[1]] ); // k = f x l
    normal = crossProduct(  M_mesh.getNodesVector()[edge->getEdge().getVerticesIds()[0]] -
                            M_mesh.getNodesVector()[edge->getEdge().getVerticesIds()[1]],
                            normal); // n = l x k
    normal.normalize();

    const Real KnDotF = fabs(dotProduct(k * normal, f)); // K n * f

    return A * KnDotF / D;
}

Real FluxOperator::findDirichletAlpha (const UInt & facetId, const Edge_ID * edge) const
{
    const Point3D borderCenter = edge->getCentroid();
    const Point3D facetCenter = M_mesh.getFacetsVector()[facetId].getCentroid();

    const Real A = edge->getSize() *
                   M_mesh.getPropertiesMap().getProperties(M_mesh.getFacetsVector()[facetId].getZoneCode()).M_aperture;

    const PermPtr_Type & k =
            M_mesh.getPropertiesMap().getProperties(M_mesh.getFacetsVector()[facetId].getZoneCode()).M_permeability;

    Point3D f = borderCenter - facetCenter;
    const Real D = sqrt(dotProduct(f, f));
    assert( D != 0 );

    f /= D;

    Point3D normal = crossProduct(  f,
                            M_mesh.getNodesVector()[edge->getEdge().getVerticesIds()[0]] -
                            M_mesh.getNodesVector()[edge->getEdge().getVerticesIds()[1]] ); // k = f x l
    normal = crossProduct(  M_mesh.getNodesVector()[edge->getEdge().getVerticesIds()[0]] -
                            M_mesh.getNodesVector()[edge->getEdge().getVerticesIds()[1]],
                            normal); // n = l x k
    normal.normalize();

    const Real KnDotF = fabs(dotProduct(k * normal, f)); // K n * f

    return A * KnDotF / D;
}

void FluxOperator::reserve_space()
{
	auto & T = *M_matrix;
	std::vector<UInt> Matrix_elements( T.cols(), 0 );
    for (auto& facet_it : M_mesh.getFractureFacetsIdsVector())
    {	
		auto f_neighbors = facet_it.getFractureNeighbors();
		UInt N_neighbors = 0;
        auto sum_maker = [&N_neighbors](std::pair< Fracture_Juncture, std::vector<UInt> > p){N_neighbors += p.second.size() ;};
        std::for_each(f_neighbors.begin(), f_neighbors.end(), sum_maker);
		Matrix_elements[ facet_it.getFractureId() ] = N_neighbors + 1;
	}	
	T.reserve(Matrix_elements);
}

void FluxOperator::assemble()
{
	auto & T = *M_matrix;
    for (auto& facet_it : M_mesh.getFractureFacetsIdsVector())
    {	
		const UInt Id_F = facet_it.getFractureId();	
		// Assemble -T matrix taking into account fractures intersections with the "star-delta" transformation
        for (auto juncture_it : facet_it.getFractureNeighbors())
        {	
            std::vector<Real> alphas;
            const Real alphaF = findFracturesAlpha (juncture_it.first, facet_it.getFractureId());

            for (auto neighbors_it : juncture_it.second)
                alphas.emplace_back(findFracturesAlpha (juncture_it.first, neighbors_it));

            Real a_sum = alphaF;
            auto sum_maker = [&a_sum](Real item){a_sum += item;};
            std::for_each(alphas.begin(), alphas.end(), sum_maker);

            for (UInt counter = 0; counter < alphas.size(); ++counter)
            {
				
                const Real QFf = alphaF * alphas[counter] * M_mesh.getPropertiesMap().getMobility() / a_sum;
                
                T.coeffRef(Id_F, Id_F) -= QFf;
                
                T.insert(Id_F,juncture_it.second[counter]) = QFf; 
			}
		}
    }
}

void FluxOperator::ImposeBConFractures(SpMat & S, Vector & rhs)
{
	auto & facetVectorRef = M_mesh.getFacetsVector();
	// assemble BC on fractures
    for (auto& edge_it : M_mesh.getBorderTipEdgesIdsVector())
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
                    borderId = (border_it > borderId) ? border_it : borderId;
            }
            else if(!isD && M_bc.getBordersBCMap().at(border_it).getBCType() == Neumann)
            {	
                borderId = (border_it > borderId) ? border_it : borderId;
			}
        }
		
		const UInt numFacetsTot = M_mesh.getFacetsVector().size() + M_mesh.getFractureFacetsIdsVector().size();
		const UInt numCells     = M_mesh.getCellsVector().size();
		
        // loop over the fracture facets to impose BC on fractures
        for(auto facet_it : edge_it.getSeparatedFacetsIds())
        {	 
            if(facetVectorRef[facet_it].isFracture())
            {	
                const UInt neighborId = facetVectorRef[facet_it].getFractureFacetId();
                const Real aperture = M_mesh.getPropertiesMap().getProperties( facetVectorRef[facet_it].getZoneCode() ).M_aperture;

                if(M_bc.getBordersBCMap().at(borderId).getBCType() == Neumann)
                {	
                    const Real Q1o = M_bc.getBordersBCMap().at(borderId).getBC()(edge_it.getCentroid()) * edge_it.getSize() * aperture;
                    rhs[ numFacetsTot+numCells + neighborId ] += Q1o;
                } // if
                
                else if(M_bc.getBordersBCMap().at(borderId).getBCType() == Dirichlet)
                {	
                    const Real alpha1 = findAlpha (facet_it, &edge_it);
                    const Real alpha2 = findDirichletAlpha (facet_it, &edge_it);

                    const Real T12 = 2 * alpha1*alpha2/(alpha1 + alpha2);
                    const Real Q12 = T12 * M_mesh.getPropertiesMap().getMobility();
                    const Real Q1o = Q12 * M_bc.getBordersBCMap().at(borderId).getBC()(edge_it.getCentroid());

                    S.coeffRef( numFacetsTot+numCells + neighborId, numFacetsTot+numCells + neighborId) -= Q12; 
                    rhs[ numFacetsTot+numCells + neighborId ] -= Q1o;
                } // else if
            } // if
        }// for
    } // for
}

void FluxOperator::ImposeBConFractures(Vector & rhs)
{
	auto & T = *M_matrix;
	auto & facetVectorRef = M_mesh.getFacetsVector();
	// assemble BC on fractures
    for (auto& edge_it : M_mesh.getBorderTipEdgesIdsVector())
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
                    borderId = (border_it > borderId) ? border_it : borderId;
            }
            else if(!isD && M_bc.getBordersBCMap().at(border_it).getBCType() == Neumann)
            {	
                borderId = (border_it > borderId) ? border_it : borderId;
			}
        }
		
		const UInt numCells     = M_mesh.getCellsVector().size();
		
        // loop over the fracture facets to impose BC on fractures
        for(auto facet_it : edge_it.getSeparatedFacetsIds())
        {	 
            if(facetVectorRef[facet_it].isFracture())
            {	
                const UInt neighborId = facetVectorRef[facet_it].getFractureFacetId();
                const Real aperture = M_mesh.getPropertiesMap().getProperties( facetVectorRef[facet_it].getZoneCode() ).M_aperture;

                if(M_bc.getBordersBCMap().at(borderId).getBCType() == Neumann)
                {	
                    const Real Q1o = M_bc.getBordersBCMap().at(borderId).getBC()(edge_it.getCentroid()) * edge_it.getSize() * aperture;
                    rhs[ numCells + neighborId ] += Q1o;
                } // if
                
                else if(M_bc.getBordersBCMap().at(borderId).getBCType() == Dirichlet)
                {	
                    const Real alpha1 = findAlpha (facet_it, &edge_it);
                    const Real alpha2 = findDirichletAlpha (facet_it, &edge_it);

                    const Real T12 = 2 * alpha1*alpha2/(alpha1 + alpha2);
                    const Real Q12 = T12 * M_mesh.getPropertiesMap().getMobility();
                    const Real Q1o = Q12 * M_bc.getBordersBCMap().at(borderId).getBC()(edge_it.getCentroid());

                    T.coeffRef( neighborId, neighborId) -= Q12; 
                    rhs[ numCells + neighborId ] -= Q1o;
                } // else if
            } // if
        }// for
    } // for
}


void global_BulkBuilder::reserve_space_Monolithic() 
{
	std::vector<UInt> Matrix_elements(S.cols(),0);
	auto & facetVectorRef   = M_mesh.getFacetsVector();
	const UInt numFacetsTot = M.cols();
	
	for (auto& cell : M_mesh.getCellsVector())
	{
	const std::vector<UInt> & cellFacetsId( cell.getFacetsIds() );
	const UInt cellId = cell.getId();
	const UInt numCellFacets = cellFacetsId.size();

	for(UInt localFacetId=0; localFacetId<numCellFacets; ++localFacetId)
    {		
		UInt globalFacetId = cellFacetsId[localFacetId];
		const Rigid_Mesh::Facet & fac = facetVectorRef[globalFacetId];
			
		// This is to take into account the decoupling of fractures facets
		if( fac.isFracture() && (cell.orientationFacet(fac)<0) )
			globalFacetId = facetVectorRef.size() + fac.getFractureFacetId();

		Matrix_elements[globalFacetId] += numCellFacets + 1;            // space for the M and B elements
		
		if( fac.isBorderFacet() &&
		( M_bc.getBordersBCMap().at(fac.getBorderId()).getBCType() == Neumann ) )
			continue;
		Matrix_elements[ numFacetsTot+cellId ] += 1;                    // space fot the Dt element
		
        }
    }
    S.reserve(Matrix_elements);	
}

void global_BulkBuilder::reserve_space() 
{
	std::vector<UInt> MMatrix_elements(M.cols(),0);
	std::vector<UInt> BMatrix_elements(B.cols(),0);
	std::vector<UInt> DtMatrix_elements(Dt.cols(),0);
	
	auto & facetVectorRef   = M_mesh.getFacetsVector();
	
	for (auto& cell : M_mesh.getCellsVector())
	{
	const std::vector<UInt> & cellFacetsId( cell.getFacetsIds() );
	const UInt cellId = cell.getId();
	const UInt numCellFacets = cellFacetsId.size();

	for(UInt localFacetId=0; localFacetId<numCellFacets; ++localFacetId)
    {		
		UInt globalFacetId = cellFacetsId[localFacetId];
		const Rigid_Mesh::Facet & fac = facetVectorRef[globalFacetId];
			
		// This is to take into account the decoupling of fractures facets
		if( fac.isFracture() && (cell.orientationFacet(fac)<0) )
			globalFacetId = facetVectorRef.size() + fac.getFractureFacetId();

		MMatrix_elements[globalFacetId] += numCellFacets;            // space for the M element
		BMatrix_elements[globalFacetId] += 1;                        // space for the B element
		
		if( fac.isBorderFacet() &&
		( M_bc.getBordersBCMap().at(fac.getBorderId()).getBCType() == Neumann ) )
			continue;
		DtMatrix_elements[cellId] += 1;    
		
        }
    }
    M.reserve(MMatrix_elements);
    B.reserve(BMatrix_elements);
    Dt.reserve(DtMatrix_elements);	
}

void global_BulkBuilder::build_Monolithic()   
{
	auto & facetVectorRef   = M_mesh.getFacetsVector();
	const UInt numFacets    = facetVectorRef.size();
	const UInt numFacetsTot = M.cols();
	
	std::cout<<"Assembling mimetic inner product and divergence..."<<std::endl;
	
	for (auto& cell : M_mesh.getCellsVector())
	{	 	
		// Definisco prodotto interno locale
		local_InnerProduct   localIP(M_mesh,cell);
		// Definisco divergenza locale
		local_Div            localDIV(M_mesh,cell);
		// Definisco builder locale
		local_builder        localBUILDER(localIP,localDIV,cell);
		// Assemblo matrici locali
		localBUILDER.build();

		auto & Mp = localIP.getMp();
		auto & Bp = localDIV.getBp();
		
		const std::vector<UInt> & cellFacetsId( cell.getFacetsIds() );
		const UInt numCellFacets  = cellFacetsId.size();
		const UInt cellId         = cell.getId();
		
		// Assemblo matrice locale in matrice globale
		for(UInt iloc=0; iloc<numCellFacets; ++iloc)
		{
            UInt i = cellFacetsId[iloc];                                 //global Id
            const Rigid_Mesh::Facet & fac = facetVectorRef[i];     
            
            // This is to take into account the decoupling of fractures facets
			if( fac.isFracture() && (cell.orientationFacet(fac)<0) )
					i = numFacets + fac.getFractureFacetId();

	        Real Bcoeff = Bp[iloc];
            if( Bcoeff != 0. )
				S.insert(numFacetsTot + cellId, i) = Bcoeff;

            Real Mcoeff = Mp(iloc,iloc);
            if( Mcoeff != 0. )
                S.coeffRef(i,i) += Mcoeff;

            for(UInt jloc=iloc+1; jloc<numCellFacets; ++jloc)
            {
                UInt j = cellFacetsId[jloc];
                const Rigid_Mesh::Facet & fac = facetVectorRef[j];
                Mcoeff = Mp(iloc,jloc);
                
                // This is to take into account the decoupling of fractures facets
				if( fac.isFracture() && (cell.orientationFacet(fac)<0) )
					j = numFacets + fac.getFractureFacetId();

                if (Mcoeff != 0.){
                    S.coeffRef(i,j) += Mcoeff;
                    S.coeffRef(j,i) += Mcoeff;
                } 
            }
            
			if( fac.isBorderFacet() &&
				( M_bc.getBordersBCMap().at(fac.getBorderId()).getBCType() == Neumann ) )
				continue;
			S.insert(i, numFacetsTot + cellId) = Bcoeff;
			
        }
	}
	std::cout<<"Done."<<std::endl;
}

void global_BulkBuilder::build()   
{
	auto & facetVectorRef   = M_mesh.getFacetsVector();
	const UInt numFacets    = facetVectorRef.size();
	
	std::cout<<"Assembling mimetic inner product and divergence..."<<std::endl;
	
	for (auto& cell : M_mesh.getCellsVector())
	{	 	
		// Definisco prodotto interno locale
		local_InnerProduct   localIP(M_mesh,cell);
		// Definisco divergenza locale
		local_Div            localDIV(M_mesh,cell);
		// Definisco builder locale
		local_builder        localBUILDER(localIP,localDIV,cell);
		// Assemblo matrici locali
		localBUILDER.build();
		
		auto & Mp = localIP.getMp();
		auto & Bp = localDIV.getBp();
		
		const std::vector<UInt> & cellFacetsId( cell.getFacetsIds() );
		const UInt numCellFacets  = cellFacetsId.size();
		const UInt cellId         = cell.getId();
		
		// Assemblo matrice locale in matrice globale
		for(UInt iloc=0; iloc<numCellFacets; ++iloc)
		{	
            UInt i = cellFacetsId[iloc];                                 //global Id
            const Rigid_Mesh::Facet & fac = facetVectorRef[i];     
            
            // This is to take into account the decoupling of fractures facets
			if( fac.isFracture() && (cell.orientationFacet(fac)<0) )
					i = numFacets + fac.getFractureFacetId();

	        Real Bcoeff = Bp[iloc];
            if( Bcoeff != 0. )
				B.insert(cellId, i) = Bcoeff;

            Real Mcoeff = Mp(iloc,iloc);
            if( Mcoeff != 0. )
                M.coeffRef(i,i) += Mcoeff;

            for(UInt jloc=iloc+1; jloc<numCellFacets; ++jloc){

                UInt j = cellFacetsId[jloc];
                const Rigid_Mesh::Facet & fac = facetVectorRef[j];
                Mcoeff = Mp(iloc,jloc);
                
                // This is to take into account the decoupling of fractures facets
				if( fac.isFracture() && (cell.orientationFacet(fac)<0) )
					j = numFacets + fac.getFractureFacetId();

                if (Mcoeff != 0.){
                    M.coeffRef(i,j) += Mcoeff;
                    M.coeffRef(j,i) += Mcoeff;
                } 
            }
            
			if( fac.isBorderFacet() &&
				( M_bc.getBordersBCMap().at(fac.getBorderId()).getBCType() == Neumann ) )
				continue;
			Dt.insert(i, cellId) = Bcoeff;
				
		}
	}
	std::cout<<"Done."<<std::endl;
}

void FractureBuilder::reserve_space_Monolithic()
{
	const UInt numFacetsTot = M.cols();
	const UInt numCells     = M_mesh.getCellsVector().size(); 
	std::vector<UInt> Matrix_elements(S.cols(),0);
    std::fill( Matrix_elements.begin() + numFacetsTot+numCells, Matrix_elements.end(), 2);
    
    for (auto& facet_it : M_mesh.getFractureFacetsIdsVector())
    {	
		//This is for -T
		auto f_neighbors = facet_it.getFractureNeighbors();
		UInt N_neighbors = 0;
        auto sum_maker = [&N_neighbors](std::pair< FluxOperator::Fracture_Juncture,std::vector<UInt> > p)
			{N_neighbors += p.second.size() ;};
        std::for_each(f_neighbors.begin(), f_neighbors.end(), sum_maker);
		Matrix_elements[ numFacetsTot+numCells + facet_it.getFractureId() ] += N_neighbors + 1;
		
		// This is for C and coupling terms on M
		UInt Id_plus   = facet_it.getId();
		UInt Id_minus  = M_mesh.getFacetsVector().size() + facet_it.getFractureId();
		Matrix_elements[ Id_plus ]   += 2;
		Matrix_elements[ Id_minus ]  += 2;	
	}
    S.reserve(Matrix_elements);
}

void FractureBuilder::reserve_space()
{
	auto & T = FO.getMatrix();
	std::vector<UInt> CMatrix_elements( C.cols(), 0 );
	std::vector<UInt> TMatrix_elements( T.cols(), 0 );
	std::vector<UInt> MMatrix_elements( M.cols(), 0 );
    for (auto& facet_it : M_mesh.getFractureFacetsIdsVector())
    {	
		//This is for -T
		auto f_neighbors = facet_it.getFractureNeighbors();
		UInt N_neighbors = 0;
        auto sum_maker = [&N_neighbors](std::pair< FluxOperator::Fracture_Juncture, std::vector<UInt> > p)
			{N_neighbors += p.second.size() ;};
        std::for_each(f_neighbors.begin(), f_neighbors.end(), sum_maker);
		TMatrix_elements[ facet_it.getFractureId() ] = N_neighbors + 1;
		
		// This is for C and coupling terms on M
		UInt Id_plus   = facet_it.getId();
		UInt Id_minus  = M_mesh.getFacetsVector().size() + facet_it.getFractureId();
		CMatrix_elements[ Id_plus ]   = 1;
		CMatrix_elements[ Id_minus ]  = 1;
		MMatrix_elements[ Id_plus ]   = 1;
		MMatrix_elements[ Id_minus ]  = 1;
	}	
	C.reserve(CMatrix_elements);
	T.reserve(TMatrix_elements);
	M.reserve(MMatrix_elements);
}

void FractureBuilder::build_Monolithic()
{
    for (auto& facet_it : M_mesh.getFractureFacetsIdsVector())
    {	
        auto& F_permeability = M_mesh.getPropertiesMap().getProperties(facet_it.getZoneCode()).M_permeability;
        const Real Kfn = F_permeability->operator()(0,0);
        const Real F_aperture = M_mesh.getPropertiesMap().getProperties(facet_it.getZoneCode()).M_aperture;
		const Real eta = F_aperture / Kfn;
		const Real xsi0     =  (2. * xsi - 1.) / 4.;
		
		const UInt numFacetsTot = M.cols();
		const UInt numCells     = M_mesh.getCellsVector().size();
		
		const UInt Id_plus         = facet_it.getId();
		const UInt Id_F            = facet_it.getFractureId();
		const UInt Id_minus        = M_mesh.getFacetsVector().size() + Id_F;
		const Real facetMeasure    = facet_it.getSize();
		
		// Coupling terms on M
		S.coeffRef(Id_plus,Id_plus)    +=   eta*facetMeasure/4. + eta*xsi0*facetMeasure;
		S.insert(Id_plus,Id_minus)      =   eta*facetMeasure/4. - eta*xsi0*facetMeasure;
		S.coeffRef(Id_minus,Id_minus)  +=   eta*facetMeasure/4. + eta*xsi0*facetMeasure;
		S.insert(Id_minus,Id_plus)      =   eta*facetMeasure/4. - eta*xsi0*facetMeasure;
		
		// Coupling matrix Ct and C
		S.insert(Id_plus, numFacetsTot+numCells + Id_F)     =  facetMeasure;
		S.insert(Id_minus, numFacetsTot+numCells + Id_F)    = -facetMeasure;
		S.insert(numFacetsTot+numCells + Id_F, Id_plus)     =  facetMeasure;
		S.insert(numFacetsTot+numCells + Id_F, Id_minus)    = -facetMeasure;
	
		// Assemble -T matrix taking into account fractures intersections with the "star-delta" transformation
        for (auto juncture_it : facet_it.getFractureNeighbors())
        {	
            std::vector<Real> alphas;
            const Real alphaF = FO.findFracturesAlpha (juncture_it.first, facet_it.getFractureId());

            for (auto neighbors_it : juncture_it.second)
                alphas.emplace_back(FO.findFracturesAlpha (juncture_it.first, neighbors_it));

            Real a_sum = alphaF;
            auto sum_maker = [&a_sum](Real item){a_sum += item;};
            std::for_each(alphas.begin(), alphas.end(), sum_maker);

            for (UInt counter = 0; counter < alphas.size(); ++counter)
            {	
                const Real QFf = alphaF * alphas[counter] * M_mesh.getPropertiesMap().getMobility() / a_sum;
                
                S.coeffRef( numFacetsTot+numCells + Id_F, numFacetsTot+numCells + Id_F ) -= QFf;
                
                S.insert( numFacetsTot+numCells + Id_F, numFacetsTot+numCells + juncture_it.second[counter]) = QFf; 
                
            }
        }
    }
}

void FractureBuilder::build()
{
	auto & T = FO.getMatrix();
    for (auto& facet_it : M_mesh.getFractureFacetsIdsVector())
    {	
        auto& F_permeability = M_mesh.getPropertiesMap().getProperties(facet_it.getZoneCode()).M_permeability;
        const Real Kfn = F_permeability->operator()(0,0);
        const Real F_aperture = M_mesh.getPropertiesMap().getProperties(facet_it.getZoneCode()).M_aperture;
		const Real eta = F_aperture / Kfn;
		const Real xsi0     =  (2. * xsi - 1.) / 4.;
		
		const UInt numFacetsTot = M.cols();
		const UInt numCells     = M_mesh.getCellsVector().size();
		
		const UInt Id_plus         = facet_it.getId();
		const UInt Id_F            = facet_it.getFractureId();
		const UInt Id_minus        = M_mesh.getFacetsVector().size() + Id_F;
		const Real facetMeasure    = facet_it.getSize();
		
		// Coupling terms on M
		M.coeffRef(Id_plus,Id_plus)    +=   eta*facetMeasure/4. + eta*xsi0*facetMeasure;
		M.insert(Id_plus,Id_minus)      =   eta*facetMeasure/4. - eta*xsi0*facetMeasure;
		M.coeffRef(Id_minus,Id_minus)  +=   eta*facetMeasure/4. + eta*xsi0*facetMeasure;
		M.insert(Id_minus,Id_plus)      =   eta*facetMeasure/4. - eta*xsi0*facetMeasure;
		
		// Coupling matrix Ct and C
		C.insert(Id_F, Id_plus)     =  facetMeasure;
		C.insert(Id_F, Id_minus)    = -facetMeasure;
	
		// Assemble -T matrix taking into account fractures intersections with the "star-delta" transformation
        for (auto juncture_it : facet_it.getFractureNeighbors())
        {	
            std::vector<Real> alphas;
            const Real alphaF = FO.findFracturesAlpha (juncture_it.first, facet_it.getFractureId());

            for (auto neighbors_it : juncture_it.second)
                alphas.emplace_back(FO.findFracturesAlpha (juncture_it.first, neighbors_it));

            Real a_sum = alphaF;
            auto sum_maker = [&a_sum](Real item){a_sum += item;};
            std::for_each(alphas.begin(), alphas.end(), sum_maker);

            for (UInt counter = 0; counter < alphas.size(); ++counter)
            {	
                const Real QFf = alphaF * alphas[counter] * M_mesh.getPropertiesMap().getMobility() / a_sum;
                
                T.coeffRef( Id_F, Id_F ) -= QFf;
                
                T.insert( Id_F, juncture_it.second[counter] ) = QFf; 
                
            }
        }
    }
}


}  //FVCode3d

