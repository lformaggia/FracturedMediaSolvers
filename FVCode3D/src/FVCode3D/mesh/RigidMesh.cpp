/*!
 * @file Rigid_Mesh.cpp
 * @brief Class for unstructured mesh (definitions).
 */

#include <FVCode3D/mesh/RigidMesh.hpp>
#include <FVCode3D/property/Properties.hpp>

namespace FVCode3D
{

// --------------------   Class Rigid_Mesh   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Rigid_Mesh (Mesh3D & generic_mesh, const PropertiesMap & prop, const bool renumber, const bool buildEdges):
    M_nodes(generic_mesh.getNodesVector()), M_properties(prop), M_renumber(renumber)
{
    std::map<UInt, UInt> oldToNewMapCells;
    std::map<UInt, UInt> oldToNewMapFacets;
    cellsVectorBuilder(generic_mesh, oldToNewMapCells);
    if(M_renumber)
        adjustCellNeighbors(oldToNewMapCells);
    facetsVectorsBuilder(generic_mesh, oldToNewMapCells, oldToNewMapFacets);
    if(M_renumber)
        adjustCellFacets(oldToNewMapFacets);
    if(buildEdges)
        edgesVectorBuilder();
}

Rigid_Mesh::Rigid_Mesh(const Rigid_Mesh & mymesh):
    M_nodes(mymesh.getNodesVector()), M_properties(mymesh.getPropertiesMap()), M_renumber(mymesh.isRenumbered())
{
    // Geometrical entities
    M_cells.reserve(mymesh.getCellsVector().size());
    for (auto& it : mymesh.getCellsVector())
    {
        M_cells.emplace_back(it, this); // Cell
    }

    M_facets.reserve(mymesh.getFacetsVector().size());
    for (auto& it : mymesh.getFacetsVector())
    {
        M_facets.emplace_back(it, this); // Facet
    }

    M_edges.reserve(mymesh.getEdgesVector().size());
    for (auto& it : mymesh.getEdgesVector())
    {
        M_edges.emplace_back(it, this); // Edge
    }

    // Facet ids vectors
    M_borderFacets.reserve(mymesh.getBorderFacetsIdsVector().size());
    for (auto& it : mymesh.getBorderFacetsIdsVector())
    {
        M_borderFacets.emplace_back(it, this); // Border_Facet
    }

    M_internalFacets.reserve(mymesh.getInternalFacetsIdsVector().size());
    for (auto& it : mymesh.getInternalFacetsIdsVector())
    {
        M_internalFacets.emplace_back(it, this); // Regular_Facet
    }

    M_fractureFacets.reserve(mymesh.getFractureFacetsIdsVector().size());
    for (auto& it : mymesh.getFractureFacetsIdsVector())
    {
        M_fractureFacets.emplace_back(it, this); // Fracture_Facet
    }

    // Edge ids vectors
    M_internalEdges.reserve(mymesh.getInternalEdgesIdsVector().size());
    for (auto& it : mymesh.getInternalEdgesIdsVector())
    {
        M_internalEdges.emplace_back(it, this); // Regular_Edge
    }

    M_pureBorderEdges.reserve(mymesh.getPureBorderEdgesIdsVector().size());
    for (auto& it : mymesh.getPureBorderEdgesIdsVector())
    {
        M_pureBorderEdges.emplace_back(it, this); // Pure_Border_Edge
    }

    M_junctureEdges.reserve(mymesh.getJunctureEdgesIdsVector().size());
    for (auto& it : mymesh.getJunctureEdgesIdsVector())
    {
        M_junctureEdges.emplace_back(it, this); // Juncture_Edge
    }

    M_internalTipEdges.reserve(mymesh.getInternalTipEdgesIdsVector().size());
    for (auto& it : mymesh.getInternalTipEdgesIdsVector())
    {
        M_internalTipEdges.emplace_back(it, this); // Internal_Tip_Edge
    }

    M_borderTipEdges.reserve(mymesh.getBorderTipEdgesIdsVector().size());
    for (auto& it : mymesh.getBorderTipEdgesIdsVector())
    {
        M_borderTipEdges.emplace_back(it, this); // Border_Tip_Edge
    }

    // Edge pointer ids vectors
    //border
    M_borderEdges.reserve(mymesh.getPureBorderEdgesIdsVector().size() + mymesh.getBorderTipEdgesIdsVector().size());
    for (auto& it : mymesh.getPureBorderEdgesIdsVector())
    {
        M_borderEdges.push_back(&it);
    }

    for (auto& it : mymesh.getBorderTipEdgesIdsVector())
    {
        M_borderEdges.push_back(&it);
    }

    //tip
    M_tipEdges.reserve(mymesh.getInternalTipEdgesIdsVector().size() + mymesh.getBorderTipEdgesIdsVector().size());
    for (auto& it : mymesh.getInternalTipEdgesIdsVector())
    {
        M_tipEdges.push_back(&it);
    }

    for (auto& it : mymesh.getBorderTipEdgesIdsVector())
    {
        M_tipEdges.push_back(&it);
    }

    //fracture
    M_fractureEdges.reserve(mymesh.getJunctureEdgesIdsVector().size() + mymesh.getTipEdgesIdsVector().size());
    for (auto& it : mymesh.getJunctureEdgesIdsVector())
    {
        M_fractureEdges.push_back(&it);
    }

    for (auto it : mymesh.getTipEdgesIdsVector())
    {
        M_fractureEdges.push_back(it);
    }
}

// ==================================================
// Protected Method
// ==================================================

void Rigid_Mesh::cellsVectorBuilder (Mesh3D & generic_mesh, std::map<UInt, UInt> & oldToNewMapCells)
{
    M_cells.reserve(generic_mesh.getCellsMap().size());

    if(!M_renumber)
    {
        for (auto& it : generic_mesh.getCellsMap())
        {
            M_cells.emplace_back(it.second, this, it.first); // Cell
        }
    }
    else
    {
        UInt position = 0;
        for (auto& it : generic_mesh.getCellsMap())
        {
            oldToNewMapCells[it.first] = position;
            M_cells.emplace_back(it.second, this, position); // Cell
            ++position;
        }
    }
}

void Rigid_Mesh::adjustCellNeighbors ( const std::map<UInt, UInt> & oldToNewMapCells )
{
    for(auto& it : M_cells)
        for(auto& Neighbors_it : it.M_neighborsIds)
            Neighbors_it = oldToNewMapCells.at(Neighbors_it);
}

void Rigid_Mesh::facetsVectorsBuilder ( Mesh3D & generic_mesh, const std::map<UInt, UInt> & oldToNewMapCells,
    std::map<UInt, UInt> & oldToNewMapFacets)
{
    // Define which facets are fracture ones, which are border ones and which are regular ones.
    const Facet3DMap & facet_container = generic_mesh.getFacetsMap();
    UInt fractureId = 0;

    M_facets.reserve(facet_container.size());
    M_fractureFacets.reserve(facet_container.size());
    M_borderFacets.reserve(facet_container.size());
    M_internalFacets.reserve(facet_container.size());

    if(!M_renumber)
    {
        for(auto& it : facet_container)
        {
            M_facets.emplace_back(it.second, this, oldToNewMapCells, it.first); // Facet

            if (it.second.isFracture())
            {
                M_fractureFacets.emplace_back(it.first, this); // Fracture_Facet

                M_facets.rbegin()->M_fractureFacetId = fractureId;
                ++fractureId;
            }
            else
            {
                if (it.second.getSeparatedCells().size() == 1)
                    M_borderFacets.emplace_back(it.first, this); // Border_Facet
                else
                    M_internalFacets.emplace_back(it.first, this); // Regular_Facet
            }
        }
    }
    else
    {
        UInt position = 0;

        for(auto& it : facet_container)
        {
            M_facets.emplace_back(it.second, this, oldToNewMapCells, position); // Facet

            oldToNewMapFacets[it.first] = position;

            if (it.second.isFracture())
            {
                M_fractureFacets.emplace_back(position, this); // Fracture_Facet

                M_facets.rbegin()->M_fractureFacetId = fractureId;
                ++fractureId;
            }
            else
            {
                if (it.second.getSeparatedCells().size() == 1)
                    M_borderFacets.emplace_back(position, this); // Border_Facet
                else
                    M_internalFacets.emplace_back(position, this); // Regular_Facet
            }
            ++ position;
        }
    }
    // TODO shrink to fit or not shrink to fit?

    // Create Junctures
    std::map<Fracture_Juncture, std::vector<UInt> > nodes_fracture_map;

    for(auto& it : M_fractureFacets)
    {
        auto vertex_it  = it.getFacet().getVerticesIds().begin();
        auto vertex_it2 = vertex_it;
        ++vertex_it2;

        for( ; vertex_it2!=it.getFacet().getVerticesIds().end(); ++vertex_it, ++vertex_it2)
        {
            //The if state is necessary in order to order the vertexes!!
            if (*(vertex_it2) < *(vertex_it))
                nodes_fracture_map[Fracture_Juncture(*(vertex_it2), *(vertex_it))].push_back(it.getFractureId());
            else
                nodes_fracture_map[Fracture_Juncture(*(vertex_it), *(vertex_it2))].push_back(it.getFractureId());
        }

        vertex_it2 = it.getFacet().getVerticesIds().begin();

        //The if state is necessary in order to order the vertexes!!
        if (*(vertex_it2) < *(vertex_it))
            nodes_fracture_map[Fracture_Juncture(*(vertex_it2), *(vertex_it))].push_back(it.getFractureId());
        else
            nodes_fracture_map[Fracture_Juncture(*(vertex_it), *(vertex_it2))].push_back(it.getFractureId());
    }

    Fracture_Juncture m_junct;

    for(auto& it : M_fractureFacets)
    {
        auto vertex_it  = it.getFacet().getVerticesIds().begin();
        auto vertex_it2 = vertex_it;
        ++vertex_it2;

        for( ; vertex_it2!=it.getFacet().getVerticesIds().end(); ++vertex_it, ++vertex_it2)
        {
            if ( *vertex_it < *vertex_it2 )
                m_junct = std::make_pair(*(vertex_it), *(vertex_it2));
            else
                m_junct = std::make_pair(*(vertex_it2), *(vertex_it));
            if (nodes_fracture_map.at(m_junct).size() > 1)
            {
                it.M_fractureNeighbors[m_junct] = nodes_fracture_map.at(m_junct);
                auto find_it = std::find(it.M_fractureNeighbors[m_junct].begin(), it.M_fractureNeighbors[m_junct].end(), it.getFractureId());
                it.M_fractureNeighbors[m_junct].erase(find_it);
            }
            else
                it.M_fractureTips.insert(static_cast<Fracture_Tip>(m_junct));
        }

        vertex_it2 = it.getFacet().getVerticesIds().begin();

        //The if state is necessary in order to order the vertexes!!
        if (*(vertex_it) < *(vertex_it2))
            m_junct = std::make_pair(*(vertex_it), *(vertex_it2));
        else
            m_junct = std::make_pair(*(vertex_it2), *(vertex_it));

        if (nodes_fracture_map.at(m_junct).size() > 1)
        {
            it.M_fractureNeighbors[m_junct] = nodes_fracture_map.at(m_junct);
            auto find_it = std::find(it.M_fractureNeighbors[m_junct].begin(), it.M_fractureNeighbors[m_junct].end(), it.getFractureId());
            it.M_fractureNeighbors[m_junct].erase(find_it);
        }
        else
            it.M_fractureTips.insert(static_cast<Fracture_Tip>(m_junct));
    }
}

void Rigid_Mesh::adjustCellFacets(const std::map<UInt, UInt> & oldToNewMapFacets)
{
    for(auto& it : M_cells)
        for(auto& Facets_it : it.M_facetsIds)
            Facets_it = oldToNewMapFacets.at(Facets_it);
}

void Rigid_Mesh::edgesBuilder()
{
    UInt edgeId = 0;
    UInt nFracFacets;

    // edge as pair -> separated facets
    std::map<Edge3D, std::vector<UInt> > edgesMap;

    // loop only over the facets
    for(auto& it : M_facets)
    {
        auto vertex_it  = it.getVerticesIds().begin();
        auto vertex_it2 = vertex_it;
        ++vertex_it2;

        for( ; vertex_it2 != it.getVerticesIds().end(); ++vertex_it, ++vertex_it2)
        {
            if (*(vertex_it) < *(vertex_it2))
            {
                edgesMap.insert(std::pair<Edge3D, std::vector<UInt> >(std::make_pair(*vertex_it, *vertex_it2), std::vector<UInt>()) );
                edgesMap[std::make_pair(*vertex_it, *vertex_it2)].push_back(it.getId());
            }
            else
            {
                edgesMap.insert(std::pair<Edge3D, std::vector<UInt> >(std::make_pair(*vertex_it2, *vertex_it), std::vector<UInt>()) );
                edgesMap[std::make_pair(*vertex_it2, *vertex_it)].push_back(it.getId());
            }
        }

        vertex_it2 = it.getVerticesIds().begin();

        if (*(vertex_it) < *(vertex_it2))
        {
            edgesMap.insert(std::pair<Edge3D, std::vector<UInt> >(std::make_pair(*vertex_it, *vertex_it2), std::vector<UInt>()) );
            edgesMap[std::make_pair(*vertex_it, *vertex_it2)].push_back(it.getId());
        }
        else
        {
            edgesMap.insert(std::pair<Edge3D, std::vector<UInt> >(std::make_pair(*vertex_it2, *vertex_it), std::vector<UInt>()) );
            edgesMap[std::make_pair(*vertex_it2, *vertex_it)].push_back(it.getId());
        }
    }

    M_edges.reserve(edgesMap.size());
    M_internalEdges.reserve(edgesMap.size());
    M_borderEdges.reserve(edgesMap.size());
    M_pureBorderEdges.reserve(edgesMap.size());
    M_fractureEdges.reserve(edgesMap.size());
    M_junctureEdges.reserve(edgesMap.size());
    M_tipEdges.reserve(edgesMap.size());
    M_internalTipEdges.reserve(edgesMap.size());
    M_borderTipEdges.reserve(edgesMap.size());

    for(auto& map_it : edgesMap)
    {
        M_edges.emplace_back(map_it.first, this, edgeId); // Edge
        M_edges.rbegin()->M_separatedFacetsIds = map_it.second;
        edgeId++;
    }

    for(auto& edge_it : M_edges)
    {
        nFracFacets = 0;
        edgeId = edge_it.getId();

        for(auto it : M_edges[edgeId].M_separatedFacetsIds)
        {
            if (M_facets[it].M_isFracture)
            {
                edge_it.M_isFracture = true;
                nFracFacets++;
            }
            edge_it.M_isBorderEdge |= M_facets[it].isBorderFacet();
        }

        edge_it.M_isTip = (nFracFacets == 1) ? true : false;
    }
}

void Rigid_Mesh::edgesIdBuilder()
{
    for(auto& edge_it : M_edges)
    {
        UInt edgeId = edge_it.getId();

        if(!edge_it.M_isFracture && !edge_it.M_isBorderEdge)
            M_internalEdges.emplace_back(edgeId,this); // Regular_Edge
        else if(!edge_it.M_isFracture && edge_it.M_isBorderEdge)
        {
            M_pureBorderEdges.emplace_back(edgeId,this); // Pure_Border_Edge
            M_borderEdges.push_back(&(*M_pureBorderEdges.rbegin()));
        }
        else if(edge_it.M_isFracture && !edge_it.M_isTip)
        {
            M_junctureEdges.emplace_back(edgeId,this); // Juncture_Edge
            M_fractureEdges.push_back(&(*M_junctureEdges.rbegin()));
        }
        else if(edge_it.M_isFracture && edge_it.M_isTip && !edge_it.M_isBorderEdge)
        {
            M_internalTipEdges.emplace_back(edgeId,this); // Internal_Tip_Edge
            M_tipEdges.push_back(&(*M_internalTipEdges.rbegin()));
            M_fractureEdges.push_back(&(*M_internalTipEdges.rbegin()));
        }
        else
        {
            M_borderTipEdges.emplace_back(edgeId,this); // Border_Tip_Edge
            M_tipEdges.push_back(&(*M_borderTipEdges.rbegin()));
            M_fractureEdges.push_back(&(*M_borderTipEdges.rbegin()));
            M_borderEdges.push_back(&(*M_borderTipEdges.rbegin()));
        }
    }
    // TODO shrink to fit or not shrink to fit?
}

bool Rigid_Mesh::hasNeighborsThroughFacet(const UInt & facet_Id, const UInt & idNeighbor) const
{
    if( M_facets[facet_Id].getSeparatedCellsIds()[0] == idNeighbor ||
       (M_facets[facet_Id].getSeparatedCellsIds().size() > 1 && M_facets[facet_Id].getSeparatedCellsIds()[1] == idNeighbor) )
    {
        return true;
    }
    return false;
}

void Rigid_Mesh::showPoint(const Point3D & generic_point, std::ostream & out) const
{
    out << generic_point.x() << " " << generic_point.y() << " " << generic_point.z() << std::endl;
}

// ==================================================
// Methods
// ==================================================

void Rigid_Mesh::showMe ( std::ostream & out ) const
{
    out << "Type = Mesh :"<< std::endl;
    out << "code = " << this << std::endl;
    out << "Number of Nodes: " << M_nodes.size() <<std::endl;
    out << "Number of Cells: " << M_cells.size() <<std::endl;
    out << "Number of Facets: " << M_facets.size() <<std::endl;
    out << "Number of Edges: " << M_edges.size() <<std::endl;
    out << "Number of Border-Facets: " << M_borderFacets.size() <<std::endl;
    out << "Number of Fracture-Facets: " << M_fractureFacets.size() <<std::endl;
    out << "Number of Standard-Facets: " << M_internalFacets.size() <<std::endl;
    out << "Number of Internal-Edges: " << M_internalEdges.size() <<std::endl;
    out << "Number of Border-Edges: " << M_borderEdges.size() <<std::endl;
    out << "Number of Pure-Border-Edges: " << M_pureBorderEdges.size() <<std::endl;
    out << "Number of Fracture-Edges: " << M_fractureEdges.size() <<std::endl;
    out << "Number of Juncture-Edges: " << M_junctureEdges.size() <<std::endl;
    out << "Number of Tip-Edges: " << M_tipEdges.size() <<std::endl;
    out << "Number of Internal-Tip-Edges: " << M_internalTipEdges.size() <<std::endl;
    out << "Number of Border-Tip-Edges: " << M_borderTipEdges.size() <<std::endl;
}

void Rigid_Mesh::BuildEdges ()
{
	for(auto & edge : M_edges)
	{
		std::vector<UInt> sepCells_byEdge;
		for(auto & facetId : edge.getSeparatedFacetsIds())
		{
			M_facets[facetId].M_edgesIds.push_back(edge.getId());
			for(auto & cellId : M_facets[facetId].getSeparatedCellsIds() ) 
				sepCells_byEdge.push_back(cellId);
		}
		sort( sepCells_byEdge.begin(), sepCells_byEdge.end() );
		sepCells_byEdge.erase( unique( sepCells_byEdge.begin(), sepCells_byEdge.end() ), sepCells_byEdge.end() );
		for(auto & cellId : sepCells_byEdge)
			M_cells[cellId].M_edgesIds.push_back(edge.getId());
	}
	// Reordinate facet edges
	for(auto & facet : M_facets)
	{
		std::vector<UInt> NewEdgeOrder;
		auto NumVertex = facet.getVerticesIds().size();
		for(auto it_vertex = facet.getVerticesIds().begin(); it_vertex!=facet.getVerticesIds().end()-1; it_vertex++)
		{
			Edge3D edge_facet(*it_vertex,*(it_vertex+1));
			for(auto & edgeId : facet.M_edgesIds)
			{
				auto & edgeVertex = M_edges[edgeId].getEdge();
				if( (edge_facet.first == edgeVertex.first && edge_facet.second == edgeVertex.second) ||
					(edge_facet.first == edgeVertex.second && edge_facet.second == edgeVertex.first) )
					NewEdgeOrder.push_back(edgeId);
			}
		}
		Edge3D edge_facet(facet.getVerticesIds()[NumVertex-1],facet.getVerticesIds()[0]);
		for(auto & edgeId : facet.M_edgesIds)
		{
			auto & edgeVertex = M_edges[edgeId].getEdge();
			if( (edge_facet.first == edgeVertex.first && edge_facet.second == edgeVertex.second) ||
				(edge_facet.first == edgeVertex.second && edge_facet.second == edgeVertex.first) )
				NewEdgeOrder.push_back(edgeId);
		}
		facet.M_edgesIds = NewEdgeOrder;		
	}
}

void Rigid_Mesh::BuildOrientationFacets ()
{
// Loop over the cells of the mesh
for(auto & cell : M_cells)
{	
	// Build map with tangent edge vectors of the facets of the cell and
	// a map to know if a facet has been already visited
	std::map<UInt, std::vector<Point3D> > tangVecToFacet;
	std::map<UInt, bool> visitedFacets;
	for(auto & facId : cell.getFacetsIds())
	{
		visitedFacets[facId] = false;
		auto & facet = M_facets[facId];
		auto NumVertex = facet.verticesNumber(); 
		for(auto it_vertex = facet.getVerticesIds().begin(); it_vertex!=facet.getVerticesIds().end()-1; it_vertex++)
		{
			tangVecToFacet[facId].push_back( M_nodes[*(it_vertex+1)] - M_nodes[*it_vertex] );
		}
		tangVecToFacet[facId].push_back( M_nodes[facet.getVerticesIds()[0]] - M_nodes[facet.getVerticesIds()[NumVertex-1]] );	
	}
	
	// Start the routine to make the orientations of the edges of the adjacent facets of the cell coherent
	std::set<UInt> adjFacIds; 
	adjFacIds.insert(cell.getFacetsIds()[0]);
	visitedFacets[cell.getFacetsIds()[0]] = true;
	UInt count = 1;
	std::set<UInt> newadjFacs;
	while(count != cell.facetsNumber())
	{
	for(auto & facId : adjFacIds)
	{
		if(count == cell.facetsNumber())
			break;
			
		auto & edgesFac = M_facets[facId].M_edgesIds;
		for(auto & edgeId : edgesFac)
		{
			if(count == cell.facetsNumber())
				break;
			
			auto facId2 = cell.adjacentFacet(M_edges[edgeId],facId);
			if(visitedFacets[facId2] == false)
			{
				auto & edgesFac2 = M_facets[facId2].M_edgesIds;
				auto & v1 = tangVecToFacet[facId][std::distance(edgesFac.begin(),std::find(edgesFac.begin(),edgesFac.end(),edgeId))];
				auto & v2 = tangVecToFacet[facId2][std::distance(edgesFac2.begin(),std::find(edgesFac2.begin(),edgesFac2.end(),edgeId))];
				Real orientationEdge = ( (dotProduct(v1,v2) >= 0.) ? 1.0 : -1.0 );
				if(orientationEdge == 1.0)
				{
					for(auto it = tangVecToFacet[facId2].begin(); it != tangVecToFacet[facId2].end(); it++)
						*it *= -1.0;
					std::reverse(tangVecToFacet[facId2].begin(),tangVecToFacet[facId2].end());
					std::reverse(edgesFac2.begin(),edgesFac2.end());
					std::reverse(M_facets[facId2].M_verticesIds.begin()+1,M_facets[facId2].M_verticesIds.end());
				}
				count++;
				visitedFacets[facId2] = true;
				newadjFacs.insert(facId2);
		    }
		}
	}
	adjFacIds.clear();
	adjFacIds = newadjFacs;
	newadjFacs.clear();
	}
	
	// Compute the normals of the facets and the signed volume of the cell, this is performed through a routine
	// for the general case of non planar facets. This is useful for possible numerical errors occurred in the
	// generation of the polyhedral mesh, i.e. "not exaclty" coplanar vertices of a facet.
	std::map<UInt, Point3D> mapNormFac;
	Real vol = 0.;
	for(auto & facetId : cell.getFacetsIds())
	{	
	auto & M_vertexIds = M_facets[facetId].M_verticesIds;
	assert(M_vertexIds.size() >= 3);
	
    if(M_vertexIds.size() > 3)
    {
        Point3D center(0., 0., 0.);
        for(auto id : M_vertexIds)
        {
            center += M_nodes[id];
        }
        center /= M_vertexIds.size();

        auto it1 = M_vertexIds.begin();
        auto it2 = it1;
        it2++;
        Real area = 0.;
        Point3D normTr(0.,0.,0.);
        mapNormFac[facetId] = Point3D(0., 0., 0.);

        for( ; it2 != M_vertexIds.end(); ++it1, ++it2)
        {
            const Point3D B = M_nodes[*it1];
            const Point3D C = M_nodes[*it2];

            area = FVCode3D::triangleArea(center, B, C);
		normTr = FVCode3D::computeNormal(center, B, C);
            mapNormFac[facetId] += normTr * area;
            vol += 1./3.*area*dotProduct(FVCode3D::triangleCentroid(center,B,C),normTr);
        }

        const Point3D B = M_nodes[ M_vertexIds[M_vertexIds.size()-1] ];
        const Point3D C = M_nodes[ M_vertexIds[0] ];

        area = FVCode3D::triangleArea( center, B, C );
        normTr = FVCode3D::computeNormal(center, B, C);

        mapNormFac[facetId] += normTr * area;
        vol += 1./3.*area*dotProduct(FVCode3D::triangleCentroid(center,B,C),normTr);
        
        mapNormFac[facetId] /= M_facets[facetId].M_area;
        mapNormFac[facetId].normalize();
    }
    else
    {
        const Point3D A(M_nodes[M_vertexIds[0]]);
        const Point3D B(M_nodes[M_vertexIds[1]]);
        const Point3D C(M_nodes[M_vertexIds[2]]);
        Point3D normTr = FVCode3D::computeNormal(A, B, C);
        Real area = FVCode3D::triangleArea( A, B, C );
        mapNormFac[facetId] = normTr;
        vol += 1./3.*area*dotProduct(FVCode3D::triangleCentroid(A,B,C),normTr);
    }
	}
	
	// Finally compute the orientation of the facets of the cell
	for(auto & facetId : cell.getFacetsIds())
	{	
		if(vol < 0.) //&& cell.getId() != 5683)
		{
			mapNormFac[facetId] *= -1.;
			cell.alphaFacets[facetId] =
				( (dotProduct(M_facets[facetId].getUnsignedNormal(),mapNormFac[facetId])>=0.) ? 1.0 : -1.0 );
		}
		else
		{
			cell.alphaFacets[facetId] =
				( (dotProduct(M_facets[facetId].getUnsignedNormal(),mapNormFac[facetId])>=0.) ? 1.0 : -1.0 );
		}
	}
	// Test to verify the volume
//	Real u = std::fabs(std::fabs(vol)-cell.getVolume())/cell.getVolume();
//	std::cout<<u<<std::endl;
}    // End loop over the cells of the mesh
}


const std::vector<Point3D> Rigid_Mesh::idsToPoints(const std::vector<UInt> & pointsIds) const
{
    std::vector<Point3D> points;
    for (auto iter : pointsIds)
        points.push_back(M_nodes[iter]);
    return points;
}

// --------------------   Class Edge   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Edge::Edge (const Edge3D & generic_edge, Rigid_Mesh * const mesh, const UInt id):
    M_edge(generic_edge), M_mesh(mesh), M_id(id), M_isBorderEdge(false), M_isFracture(false), M_isTip(false){}

Rigid_Mesh::Edge::Edge (const Edge & edge, Rigid_Mesh * const mesh):
    M_edge(edge.getEdge()), M_mesh(mesh), M_id(edge.getId()),
    M_separatedFacetsIds(edge.getSeparatedFacetsIds()), M_isBorderEdge(edge.isBorderEdge()),
    M_isFracture(edge.isFracture()), M_isTip(edge.isTip()){}

Rigid_Mesh::Edge::Edge (const Edge & edge):
    M_edge(edge.getEdge()), M_mesh(edge.getMesh()), M_id(edge.getId()),
    M_separatedFacetsIds(edge.getSeparatedFacetsIds()), M_isBorderEdge(edge.isBorderEdge()),
    M_isFracture(edge.isFracture()), M_isTip(edge.isTip()){}

const Point3D Rigid_Mesh::Edge::getCentroid () const
{
    Point3D first (M_mesh->getNodesVector()[M_edge.first]);
    Point3D second (M_mesh->getNodesVector()[M_edge.second]);
    return (first+second) / 2.;
}

Real Rigid_Mesh::Edge::length() const
{
    const Point3D & first (M_mesh->getNodesVector()[M_edge.first]);
    const Point3D & second (M_mesh->getNodesVector()[M_edge.second]);
    return (second-first).norm();
}

void Rigid_Mesh::Edge::showMe (std::ostream & out) const
{
    out << "Type = Edge: " << std::endl;
    out << "ID: " << M_id << std::endl;
    out << "M_mesh = " << M_mesh << std::endl;
    out << "Edge: (" << M_edge.first << "," << M_edge.second << ")" << std::endl;
    out << "Size: " << length() << std::endl;
    out << "Center: ";
    M_mesh->showPoint(getCentroid(), out);
    out << "Is border? " << M_isBorderEdge << std::endl;
    out << "Is fracture? " << M_isFracture << std::endl;
    out << "Is tip? " << M_isTip << std::endl;
    out << "M_separatedFacets : size = " << M_separatedFacetsIds.size();
    out << "    [ ";
    if( !M_separatedFacetsIds.empty() )
        for(auto it : M_separatedFacetsIds)
            out << it << " ";
    out << "] " << std::endl;
}

// --------------------   Class Facet   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Facet::Facet(const Facet3D & generic_facet, Rigid_Mesh * const mesh,
                         const std::map<UInt,UInt> & oldToNewMapCells, const UInt m_id):
    M_mesh(mesh), M_id(m_id), M_verticesIds(generic_facet.getVerticesVector()),
    M_area(generic_facet.getArea()),
    M_centroid(generic_facet.getCentroid()), M_unsignedNormal(generic_facet.getUnsignedNormal()),
    M_isFracture(generic_facet.isFracture()), M_fractureFacetId(0), M_borderId(generic_facet.getBorderId()),
    M_representedFractureIds(generic_facet.getRepresentedFractureIds()), M_zone(generic_facet.getZoneCode())
{
    if(!M_mesh->M_renumber)
        for(auto it : generic_facet.getSeparatedCells())
            M_separatedCellsIds.push_back(it);
    else
        for(auto it : generic_facet.getSeparatedCells())
            M_separatedCellsIds.push_back(oldToNewMapCells.at(it));
}

Rigid_Mesh::Facet::Facet(const Facet & facet):
    M_mesh(facet.getMesh()), M_id(facet.getId()), M_verticesIds(facet.getVerticesIds()), M_edgesIds(facet.getEdgesIds()),
    M_separatedCellsIds(facet.getSeparatedCellsIds()), M_area(facet.area()), M_centroid(facet.getCentroid()),
    M_unsignedNormal(facet.getUnsignedNormal()),
    M_isFracture(facet.isFracture()), M_fractureFacetId(facet.getFractureFacetId()), M_borderId(facet.getBorderId()),
    M_representedFractureIds(facet.getRepresentedFractureIds()), M_zone(facet.getZoneCode()){}

Rigid_Mesh::Facet::Facet(const Facet & facet, Rigid_Mesh * const mesh):
    M_mesh(mesh), M_id(facet.getId()), M_verticesIds(facet.getVerticesIds()), M_edgesIds(facet.getEdgesIds()),
    M_separatedCellsIds(facet.getSeparatedCellsIds()), M_area(facet.area()), M_centroid(facet.getCentroid()),
    M_unsignedNormal(facet.getUnsignedNormal()),
    M_isFracture(facet.isFracture()), M_fractureFacetId(facet.getFractureFacetId()), M_borderId(facet.getBorderId()),
    M_representedFractureIds(facet.getRepresentedFractureIds()), M_zone(facet.getZoneCode()){}

// ==================================================
// Methods
// ==================================================

void Rigid_Mesh::Facet::showMe(std::ostream & out) const
{
    out << "Type = Facet: " << std::endl;
    out << "ID: " << M_id << std::endl;
    out << "M_mesh = " << M_mesh << std::endl;
    out << "Number of Nodes: " << M_verticesIds.size()<<" : "<<std::endl;
    UInt counter=1;
    for(auto j : M_verticesIds)
    {
        out << "M_idP"<<counter<<" = " <<j<<"\t";
        ++counter;
    }
    out << std::endl;
    out << "Size: " << M_area << std::endl;
    out << "Center: ";
    M_mesh->showPoint(M_centroid, out);
    out << "M_separatedCells : size = " << M_separatedCellsIds.size();
    out << "    [ ";
    if( !M_separatedCellsIds.empty() )
        for(auto it : M_separatedCellsIds)
            out << it << " ";
    out << "] " << std::endl;
}

// --------------------   Class Cell   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Cell::Cell( const Cell3D & generic_cell, Rigid_Mesh * const mesh, const UInt m_id):
    M_mesh(mesh), M_Id(m_id), M_zoneCode(generic_cell.getZoneCode()),
    M_verticesIds(generic_cell.getVerticesVector()),
    M_centroid(generic_cell.getCentroid()),
    M_volume(generic_cell.volume())
{
    for(auto it : generic_cell.getFacetsSet())
        M_facetsIds.push_back(it);
    for(auto it : generic_cell.getNeighborsSet())
        M_neighborsIds.push_back(it);
}

Rigid_Mesh::Cell::Cell(const Cell & cell, Rigid_Mesh * const mesh):
    M_mesh(mesh), M_Id(cell.getId()), M_zoneCode(cell.getZoneCode()),
    M_verticesIds(cell.getVerticesIds()), M_edgesIds(cell.getEdgesIds()),
    M_facetsIds(cell.getFacetsIds()), M_neighborsIds(cell.getNeighborsIds()),
    M_centroid(cell.getCentroid()), M_volume(cell.getVolume()){}

Rigid_Mesh::Cell::Cell(const Cell & cell):
    M_mesh(cell.getMesh()), M_Id(cell.getId()), M_zoneCode(cell.getZoneCode()),
    M_verticesIds(cell.getVerticesIds()), M_edgesIds(cell.getEdgesIds()),
    M_facetsIds(cell.getFacetsIds()), M_neighborsIds(cell.getNeighborsIds()),
    M_centroid(cell.getCentroid()), M_volume(cell.getVolume()){}

Real Rigid_Mesh::Cell::getAlpha(const UInt facetId) const
{ 	
	auto FacetIsInCell = find(M_facetsIds.begin(),M_facetsIds.end(),facetId);
	if(FacetIsInCell == M_facetsIds.end())
	{
		std::stringstream error;
		error << "Error in getAlpha method: the facet is not in the cell.";
		throw std::runtime_error(error.str());
	}
	return alphaFacets.at(facetId);	
}

Real Rigid_Mesh::Cell::orientationFacet_forConvexCells(const Facet & fac) const
{
    bool found(false);
    UInt facetId = fac.getId();

    // control that facetId is a facet of this cell
    for( std::vector<UInt>::const_iterator it = M_facetsIds.begin(); it != M_facetsIds.end(); ++it )
    {
        if( *it == facetId )
            found = true;
    }

    // check if the facet exists in this cell
    if( !found )
    {
        std::stringstream error;
        error << "Error: the facet " << facetId << " is not in the cell.";
        throw std::runtime_error(error.str());
    }

    const Point3D facetNormal = fac.getUnsignedNormal();
    const Point3D g           = fac.getCentroid() - Rigid_Mesh::Cell::M_centroid;
    
	return ( (dotProduct(g,facetNormal)>=0.) ? 1.0 : -1.0 );	
}

// ==================================================
// Methods
// ==================================================

void Rigid_Mesh::Cell::showMe(std::ostream  & out) const
{
    out << "Type = Cell :" << std::endl;
    out << "ID: " << M_Id << std::endl;
    out << "M_mesh = " << M_mesh << std::endl;
    out << "Number of Nodes: " << M_verticesIds.size()<<" : "<<std::endl;
    UInt counter=1;
    for(auto j : M_verticesIds)
    {
        out << "M_idP"<<counter<<" = " <<j<<"\t";
        ++counter;
    }
    out << std::endl;
    out << "Volume: " << M_volume << std::endl;
    out << "Centroid: ";
    M_mesh->showPoint(M_centroid, out);
    out << "M_Neighbors_Ids : size = " << M_neighborsIds.size();
    out << "    [ ";
    if( !M_neighborsIds.empty() )
        for(auto it : M_neighborsIds)
            out << it << " ";
    out << "] " << std::endl;
}

UInt Rigid_Mesh::Cell::adjacentFacet (const Edge & edge, const UInt IdFac) const
{
	auto edgeIsInCell = find(M_edgesIds.begin(),M_edgesIds.end(),edge.getId());
	if(edgeIsInCell == M_edgesIds.end())
	{
		std::stringstream error;
		error << "Error in adjacentFacet method: the edge is not in the cell.";
		throw std::runtime_error(error.str());		
	}
	auto facetIsInCell = find(M_facetsIds.begin(),M_facetsIds.end(),IdFac);
	if(facetIsInCell == M_facetsIds.end())
	{
		std::stringstream error;
		error << "Error in adjacentFacet method: the facet is not in the cell.";
		throw std::runtime_error(error.str());		
	}
	auto edgeIsInFacet = find(M_mesh->getFacetsVector()[IdFac].getEdgesIds().begin(),
		M_mesh->getFacetsVector()[IdFac].getEdgesIds().end(),edge.getId());
	if(edgeIsInFacet == M_mesh->getFacetsVector()[IdFac].getEdgesIds().end())
	{
		std::stringstream error;
		error << "Error in adjacentFacet method: the edge is not in the facet.";
		throw std::runtime_error(error.str());		
	}
	UInt adjFacet;
	for(auto & facetId:edge.getSeparatedFacetsIds())
	{
		auto FacetIt = find(M_facetsIds.begin(),M_facetsIds.end(),facetId);
		if( FacetIt != M_facetsIds.end() && facetId != IdFac )
		{
			adjFacet = facetId;
			break;
		}
	}
	return adjFacet;
}

// --------------------   Class Edge_ID   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Edge_ID::Edge_ID(const UInt edge_id, Rigid_Mesh * const mesh):
    M_id(edge_id), M_mesh(mesh){}

// --------------------   Class Regular_Edge   --------------------

Rigid_Mesh::Regular_Edge::Regular_Edge (const UInt edge_Id, Rigid_Mesh * const mesh):
    Edge_ID(edge_Id, mesh) {}

Rigid_Mesh::Regular_Edge::Regular_Edge (const Regular_Edge & regular_edge, Rigid_Mesh * const mesh):
    Edge_ID(regular_edge.getId(), mesh) {}

Rigid_Mesh::Regular_Edge::Regular_Edge (const Regular_Edge & e):
    Edge_ID(e.getId(), e.getMesh()) {}

// --------------------   Class Border_Edge   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Border_Edge::Border_Edge (const UInt edge_Id, Rigid_Mesh * const mesh):
    Edge_ID(edge_Id, mesh)
{
    for(auto facet_it : mesh->getEdgesVector()[edge_Id].getSeparatedFacetsIds())
    {
        M_borderIds.insert(mesh->getFacetsVector()[facet_it].getBorderId());
    }
    M_borderIds.erase(0);
}

Rigid_Mesh::Border_Edge::Border_Edge (const Border_Edge & border_edge, Rigid_Mesh * const mesh):
    Edge_ID(border_edge.getId(), mesh), M_borderIds(border_edge.getBorderIds()){}

Rigid_Mesh::Border_Edge::Border_Edge (const Border_Edge & e):
    Edge_ID(e.getId(), e.getMesh()), M_borderIds(e.getBorderIds()){}

// --------------------   Class Pure_Border_Edge   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Pure_Border_Edge::Pure_Border_Edge (const UInt edge_Id, Rigid_Mesh * const mesh):
    Edge_ID(edge_Id, mesh),
    Border_Edge(edge_Id, mesh){}

Rigid_Mesh::Pure_Border_Edge::Pure_Border_Edge (const Pure_Border_Edge & pure_border_edge, Rigid_Mesh * const mesh):
    Edge_ID(pure_border_edge.getId(), mesh),
    Border_Edge(pure_border_edge.getId(), mesh){}

Rigid_Mesh::Pure_Border_Edge::Pure_Border_Edge (const Pure_Border_Edge & e):
    Edge_ID(e.getId(), e.getMesh()),
    Border_Edge(e.getId(), e.getMesh()){}

// --------------------   Class Fracture_Edge   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Fracture_Edge::Fracture_Edge (const UInt edge_Id, Rigid_Mesh * const mesh):
    Edge_ID(edge_Id, mesh)
{
    for(auto facet_it : mesh->getEdgesVector()[edge_Id].getSeparatedFacetsIds())
    {
        for(auto fracture_it : mesh->getFacetsVector()[facet_it].getRepresentedFractureIds())
            M_fractureIds.insert(fracture_it);
    }
}

Rigid_Mesh::Fracture_Edge::Fracture_Edge (const Fracture_Edge & fracture_edge, Rigid_Mesh * const mesh):
    Edge_ID(fracture_edge.getId(), mesh), M_fractureIds(fracture_edge.getFractureIds()){}

Rigid_Mesh::Fracture_Edge::Fracture_Edge (const Fracture_Edge & e):
    Edge_ID(e.getId(), e.getMesh()), M_fractureIds(e.getFractureIds()){}

// --------------------   Class Juncture_Edge   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Juncture_Edge::Juncture_Edge (const UInt edge_Id, Rigid_Mesh * const mesh):
    Edge_ID(edge_Id, mesh),
    Fracture_Edge(edge_Id, mesh){}

Rigid_Mesh::Juncture_Edge::Juncture_Edge (const Juncture_Edge & juncture_edge, Rigid_Mesh * const mesh):
    Edge_ID(juncture_edge.getId(), mesh),
    Fracture_Edge(juncture_edge.getId(), mesh){}

Rigid_Mesh::Juncture_Edge::Juncture_Edge (const Juncture_Edge & e):
    Edge_ID(e.getId(), e.getMesh()),
    Fracture_Edge(e.getId(), e.getMesh()){}

// --------------------   Class Tip_Edge   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Tip_Edge::Tip_Edge (const UInt edge_Id, Rigid_Mesh * const mesh):
    Edge_ID(edge_Id, mesh),
    Fracture_Edge(edge_Id, mesh){}

Rigid_Mesh::Tip_Edge::Tip_Edge (const Tip_Edge & tip_edge, Rigid_Mesh * const mesh):
    Edge_ID(tip_edge.getId(), mesh),
    Fracture_Edge(tip_edge.getId(), mesh){}

Rigid_Mesh::Tip_Edge::Tip_Edge (const Tip_Edge & e):
    Edge_ID(e.getId(), e.getMesh()),
    Fracture_Edge(e.getId(), e.getMesh()){}

// --------------------   Class Interior_Tip_Edge   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Internal_Tip_Edge::Internal_Tip_Edge (const UInt edge_Id, Rigid_Mesh * const mesh):
    Edge_ID(edge_Id, mesh),
    Fracture_Edge(edge_Id, mesh),
    Tip_Edge(edge_Id, mesh){}

Rigid_Mesh::Internal_Tip_Edge::Internal_Tip_Edge (const Internal_Tip_Edge & internal_tip_edge, Rigid_Mesh * const mesh):
    Edge_ID(internal_tip_edge.getId(), mesh),
    Fracture_Edge(internal_tip_edge.getId(), mesh),
    Tip_Edge(internal_tip_edge.getId(), mesh){}

Rigid_Mesh::Internal_Tip_Edge::Internal_Tip_Edge (const Internal_Tip_Edge & e):
    Edge_ID(e.getId(), e.getMesh()),
    Fracture_Edge(e.getId(), e.getMesh()),
    Tip_Edge(e.getId(), e.getMesh()){}

// --------------------   Class Border_Tip_Edge   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Border_Tip_Edge::Border_Tip_Edge (const UInt edge_Id, Rigid_Mesh * const mesh):
    Edge_ID(edge_Id, mesh),
    Fracture_Edge(edge_Id, mesh),
    Tip_Edge(edge_Id, mesh),
    Border_Edge(edge_Id, mesh){}

Rigid_Mesh::Border_Tip_Edge::Border_Tip_Edge (const Border_Tip_Edge & border_tip_edge, Rigid_Mesh * const mesh):
    Edge_ID(border_tip_edge.getId(), mesh),
    Fracture_Edge(border_tip_edge.getId(), mesh),
    Tip_Edge(border_tip_edge.getId(), mesh),
    Border_Edge(border_tip_edge.getId(), mesh){}

Rigid_Mesh::Border_Tip_Edge::Border_Tip_Edge (const Border_Tip_Edge & e):
    Edge_ID(e.getId(), e.getMesh()),
    Fracture_Edge(e.getId(), e.getMesh()),
    Tip_Edge(e.getId(), e.getMesh()),
    Border_Edge(e.getId(), e.getMesh()){}

// --------------------   Class Facet_ID   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Facet_ID::Facet_ID(const UInt facet_Id, Rigid_Mesh * const mesh):
    M_id(facet_Id), M_mesh(mesh){}

// --------------------   Class Regular_Facet   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Regular_Facet::Regular_Facet(const UInt facet_Id, Rigid_Mesh * const mesh):
    Facet_ID(facet_Id, mesh){}

Rigid_Mesh::Regular_Facet::Regular_Facet(const Regular_Facet & regular_facet):
    Facet_ID(regular_facet.getId(), regular_facet.getMesh()) {}

Rigid_Mesh::Regular_Facet::Regular_Facet(const Regular_Facet & regular_facet, Rigid_Mesh * const mesh):
    Facet_ID(regular_facet.getId(), mesh) {}

// --------------------   Class Border_Facet   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Border_Facet::Border_Facet(const UInt facet_Id, Rigid_Mesh * const mesh):
    Facet_ID(facet_Id, mesh){}

Rigid_Mesh::Border_Facet::Border_Facet(const Border_Facet & border_facet, Rigid_Mesh * const mesh):
    Facet_ID(border_facet.getId(), mesh){}

Rigid_Mesh::Border_Facet::Border_Facet(const Border_Facet & border_facet):
    Facet_ID(border_facet.getId(), border_facet.getMesh()) {}

// --------------------   Class Fracture_Facet   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Fracture_Facet::Fracture_Facet(const UInt facet_Id, Rigid_Mesh * const mesh):
    Facet_ID(facet_Id, mesh) {}

Rigid_Mesh::Fracture_Facet::Fracture_Facet(const Fracture_Facet & fracture_facet, Rigid_Mesh * const mesh):
    Facet_ID(fracture_facet.getId(), mesh),
    M_fractureNeighbors(fracture_facet.getFractureNeighbors()),
    M_fractureTips(fracture_facet.getFractureTips()) {}

Rigid_Mesh::Fracture_Facet::Fracture_Facet(const Fracture_Facet & fracture_facet):
    Facet_ID(fracture_facet.getId(), fracture_facet.getMesh()),
    M_fractureNeighbors(fracture_facet.getFractureNeighbors()),
    M_fractureTips(fracture_facet.getFractureTips()) {}

} // namespace FVCode3D
