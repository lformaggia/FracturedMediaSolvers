 /*!
 * @file Rigid_Mesh.cpp
 * @brief Class for unstructured mesh (definitions).
 */

#include "mesh/Rigid_Mesh.hpp"
#include "mesh/Properties.hpp"

namespace Geometry{

// --------------------   Class Rigid_Mesh   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Rigid_Mesh (Generic_Mesh & generic_mesh, const PropertiesMap & prop, const bool renumber):
	M_nodes(generic_mesh.getNodesVector()), M_properties(prop), M_renumber(renumber)
{
	std::map<UInt, UInt> old_to_new_mapCells;
	std::map<UInt, UInt> old_to_new_mapFacets;
	CellsVectorBuilder(generic_mesh, old_to_new_mapCells);
	if(M_renumber)
		AdjustCellNeighbors(old_to_new_mapCells);
	FacetsVectorsBuilder(generic_mesh, old_to_new_mapCells, old_to_new_mapFacets);
	if(M_renumber)
		AdjustCellFacets(old_to_new_mapFacets);
	EdgesVectorBuilder();
}

Rigid_Mesh::Rigid_Mesh(const Rigid_Mesh & mymesh):
	M_nodes(mymesh.getNodesVector()), M_properties(mymesh.getPropertiesMap()), M_renumber(mymesh.isRenumbered())
{
	for (auto it = mymesh.getCellsVector().begin(); it != mymesh.getCellsVector().end(); ++it)
	{
		Cell my_cell (*(it), this);
		M_cells.push_back (my_cell);
	}

	for (auto it = mymesh.getFacetsVector().begin(); it != mymesh.getFacetsVector().end(); ++it)
	{
		Facet my_facet (*(it), this);
		M_facets.push_back (my_facet);
	}

	for (auto it = mymesh.getBorderFacetsIdsVector().begin(); it != mymesh.getBorderFacetsIdsVector().end(); ++it)
	{
		M_borderFacets.emplace_back (Border_Facet (*(it), this));
	}

	for (auto it = mymesh.getInternalFacetsIdsVector().begin(); it != mymesh.getInternalFacetsIdsVector().end(); ++it)
	{
		M_internalFacets.emplace_back (Regular_Facet (*(it), this));
	}

	for (auto it = mymesh.getFractureFacetsIdsVector().begin(); it != mymesh.getFractureFacetsIdsVector().end(); ++it)
	{
		M_fractureFacets.emplace_back (Fracture_Facet (*(it), this));
	}
}

// ==================================================
// Protected Method
// ==================================================

void Rigid_Mesh::CellsVectorBuilder (Generic_Mesh & generic_mesh, std::map<UInt, UInt> & old_to_new_map)
{
	if(!M_renumber)
	{
		for (auto it = generic_mesh.getCellsMap().begin(); it != generic_mesh.getCellsMap().end(); ++it)
		{
			M_cells.push_back(Cell(it->second, this, it->first));
		}
	}
	else
	{
		UInt position = 0;
		for (auto it = generic_mesh.getCellsMap().begin(); it != generic_mesh.getCellsMap().end(); ++it)
		{
			old_to_new_map[it->first] = position;
			M_cells.push_back(Cell(it->second, this, position));
			++position;
		}
	}

}

void Rigid_Mesh::AdjustCellNeighbors ( const std::map<UInt, UInt> & old_to_new_map )
{
	for (auto it = M_cells.begin(); it != M_cells.end(); ++it)
		for (auto Neighbors_it = it->M_neighborsIds.begin(); Neighbors_it != it->M_neighborsIds.end(); ++ Neighbors_it)
			*(Neighbors_it) = old_to_new_map.at(*(Neighbors_it));
}

void Rigid_Mesh::FacetsVectorsBuilder ( Generic_Mesh & generic_mesh, const std::map<UInt, UInt> & old_to_new_mapCells, std::map<UInt, UInt> & old_to_new_mapFacets)
{
	// Define which facets are fracture ones, which are border ones and which are regular ones.
	M_facetsVectorsBuilder(old_to_new_mapCells, old_to_new_mapFacets, generic_mesh);

	// Create Junctures
	std::map<Fracture_Juncture, std::vector<UInt> > nodes_fracture_map;

	for(auto it = M_fractureFacets.begin(); it != M_fractureFacets.end(); ++it)
	{
		auto vertex_it  = it->getFacet().getVertexesIds().begin();
		auto vertex_it2 = vertex_it;
		++vertex_it2;

		for( ; vertex_it2!=it->getFacet().getVertexesIds().end(); ++vertex_it, ++vertex_it2)
		{
			//The if state is necessary in order to order the vertexes!!
			if (*(vertex_it2) < *(vertex_it))
				nodes_fracture_map[Fracture_Juncture(*(vertex_it2), *(vertex_it))].emplace_back(it->getId());
			else
				nodes_fracture_map[Fracture_Juncture(*(vertex_it), *(vertex_it2))].emplace_back(it->getId());
		}

		vertex_it2 = it->getFacet().getVertexesIds().begin();

		//The if state is necessary in order to order the vertexes!!
		if (*(vertex_it2) < *(vertex_it))
			nodes_fracture_map[Fracture_Juncture(*(vertex_it2), *(vertex_it))].emplace_back(it->getId());
		else
			nodes_fracture_map[Fracture_Juncture(*(vertex_it), *(vertex_it2))].emplace_back(it->getId());
	}

	Fracture_Juncture m_junct;

	for (auto it = M_fractureFacets.begin(); it != M_fractureFacets.end(); ++it)
	{
		auto vertex_it  = it->getFacet().getVertexesIds().begin();
		auto vertex_it2 = vertex_it;
		++vertex_it2;

		for( ; vertex_it2!=it->getFacet().getVertexesIds().end(); ++vertex_it, ++vertex_it2)
		{
			if ( *vertex_it < *vertex_it2 )
				m_junct = std::make_pair(*(vertex_it), *(vertex_it2));
			else
				m_junct = std::make_pair(*(vertex_it2), *(vertex_it));
			if (nodes_fracture_map.at(m_junct).size() > 1)
			{
				it->Fracture_Neighbors[m_junct] = nodes_fracture_map.at(m_junct);
				auto find_it = std::find(it->Fracture_Neighbors[m_junct].begin(), it->Fracture_Neighbors[m_junct].end(), it->getId());
				it->Fracture_Neighbors[m_junct].erase(find_it);
			}
			else
				it->Fracture_Tips.insert(static_cast<Fracture_Tip>(m_junct));
		}

		vertex_it2 = it->getFacet().getVertexesIds().begin();

		//The if state is necessary in order to order the vertexes!!
		if (*(vertex_it) < *(vertex_it2))
			m_junct = std::make_pair(*(vertex_it), *(vertex_it2));
		else
			m_junct = std::make_pair(*(vertex_it2), *(vertex_it));

		if (nodes_fracture_map.at(m_junct).size() > 1)
		{
			it->Fracture_Neighbors[m_junct] = nodes_fracture_map.at(m_junct);
			auto find_it = std::find(it->Fracture_Neighbors[m_junct].begin(), it->Fracture_Neighbors[m_junct].end(), it->getId());
			it->Fracture_Neighbors[m_junct].erase(find_it);
		}
		else
			it->Fracture_Tips.insert(static_cast<Fracture_Tip>(m_junct));
	}
}

void Rigid_Mesh::M_facetsVectorsBuilder(const std::map<UInt,UInt> & old_to_new_mapCells, std::map<UInt,UInt> & old_to_new_mapFacets, Generic_Mesh & generic_mesh)
{
	UInt fractureId = 0;

	const Generic_FacetsContainer & facet_container = generic_mesh.getFacetsMap();
	//std::map<Generic_Edge, std::vector<UInt> > edgesMap;

	if(!M_renumber)
	{
		for(auto it = facet_container.begin(); it != facet_container.end(); ++it)
		{
			M_facets.push_back(Facet( it->second, this, old_to_new_mapCells, it->first));

			if (it->second.isFracture())
			{
				std::tuple<UInt,UInt,UInt,UInt> fracture_Ids;
				std::get<0>(fracture_Ids) = it->first;
				std::get<1>(fracture_Ids) = fractureId;
				std::get<2>(fracture_Ids) = M_cells.size();
				std::get<3>(fracture_Ids) = it->second.getZoneCode();

				M_fractureFacets.push_back(
					Fracture_Facet(fracture_Ids, it->second.getRepresentedFractureIds(), this)
				);
				++fractureId;
			}
			else
			{
				if (it->second.getSeparatedCells().size() == 1)
					M_borderFacets.push_back(
						Border_Facet(it->first, it->second.getZoneCode(), it->second.getBorderId(), this)
					);
				else
					M_internalFacets.emplace_back(Regular_Facet (it->first, it->second.getZoneCode(), this));
			}
		}
	}
	else
	{
		UInt position = 0;

		for(auto it = facet_container.begin(); it != facet_container.end(); ++it)
		{
			Facet facet( it->second, this, old_to_new_mapCells, position);
			M_facets.push_back(facet);

			old_to_new_mapFacets[it->first] = position;

			if (it->second.isFracture())
			{
				std::tuple<UInt,UInt,UInt,UInt> fracture_Ids;
				std::get<0>(fracture_Ids) = position;
				std::get<1>(fracture_Ids) = fractureId;
				std::get<2>(fracture_Ids) = M_cells.size();
				std::get<3>(fracture_Ids) = it->second.getZoneCode();

				M_fractureFacets.push_back(
					Fracture_Facet( fracture_Ids, it->second.getRepresentedFractureIds(), this)
				);
/*
				auto vertex_it  = it->second.getVertexesVector().begin();
				auto vertex_it2 = vertex_it;
				++vertex_it2;

				for( ; vertex_it2 != it->second.getVertexesVector().end(); ++vertex_it, ++vertex_it2)
				{
					if (*(vertex_it) < *(vertex_it2))
					{
						edgesMap.insert(std::pair<Generic_Edge, std::vector<UInt> >(std::make_pair(*vertex_it, *vertex_it2), std::vector<UInt>()) );
						edgesMap[std::make_pair(*vertex_it, *vertex_it2)].push_back(it->getFacetId());
					}
					else
					{
						edgesMap.insert(std::pair<Generic_Edge, std::vector<UInt> >(std::make_pair(*vertex_it2, *vertex_it), std::vector<UInt>()) );
						edgesMap[std::make_pair(*vertex_it2, *vertex_it)].push_back(it->getFacetId());
					}
				}

				vertex_it2 = it->getFacet().getVertexesIds().begin();

				if (*(vertex_it) < *(vertex_it2))
				{
					edgesMap.insert(std::pair<Generic_Edge, std::vector<UInt> >(std::make_pair(*vertex_it, *vertex_it2), std::vector<UInt>()) );
					edgesMap[std::make_pair(*vertex_it, *vertex_it2)].push_back(it->getFacetId());
				}
				else
				{
					edgesMap.insert(std::pair<Generic_Edge, std::vector<UInt> >(std::make_pair(*vertex_it2, *vertex_it), std::vector<UInt>()) );
					edgesMap[std::make_pair(*vertex_it2, *vertex_it)].push_back(it->getFacetId());
				}*/

				++fractureId;
			}
			else
			{
				if (it->second.getSeparatedCells().size() == 1)
				{
					Border_Facet border_facet (position, it->second.getZoneCode(), it->second.getBorderId(), this);
					M_borderFacets.push_back (border_facet);
				}
				else
					M_internalFacets.emplace_back(Regular_Facet (position, it->second.getZoneCode(), this));
			}
			++ position;
		}
	}
}

void Rigid_Mesh::AdjustCellFacets(const std::map<UInt, UInt> & old_to_new_mapFacets)
{
	for (auto it = M_cells.begin(); it != M_cells.end(); ++it)
		for (auto Facets_it = it->M_facetsIds.begin(); Facets_it != it->M_facetsIds.end(); ++ Facets_it)
			*(Facets_it) = old_to_new_mapFacets.at(*(Facets_it));
}

void Rigid_Mesh::EdgesVectorBuilder()
{
	/*
	UInt edgeId = 0;

	std::map<Generic_Edge, std::vector<UInt> > edgesMap;
	std::set<UInt> sepFacets;

	// loop only over the border facet! Modify this if you want to store all edges
	for(auto it = M_borderFacets.begin(); it != M_borderFacets.end(); ++it)
	{
		auto vertex_it  = it->getFacet().getVertexesIds().begin();
		auto vertex_it2 = vertex_it;
		++vertex_it2;

		for( ; vertex_it2 != it->getFacet().getVertexesIds().end(); ++vertex_it, ++vertex_it2)
		{
			if (*(vertex_it) < *(vertex_it2))
			{
				edgesMap.insert(std::pair<Generic_Edge, std::vector<UInt> >(std::make_pair(*vertex_it, *vertex_it2), std::vector<UInt>()) );
				edgesMap[std::make_pair(*vertex_it, *vertex_it2)].push_back(it->getFacetId());
			}
			else
			{
				edgesMap.insert(std::pair<Generic_Edge, std::vector<UInt> >(std::make_pair(*vertex_it2, *vertex_it), std::vector<UInt>()) );
				edgesMap[std::make_pair(*vertex_it2, *vertex_it)].push_back(it->getFacetId());
			}
		}

		vertex_it2 = it->getFacet().getVertexesIds().begin();

		if (*(vertex_it) < *(vertex_it2))
		{
			edgesMap.insert(std::pair<Generic_Edge, std::vector<UInt> >(std::make_pair(*vertex_it, *vertex_it2), std::vector<UInt>()) );
			edgesMap[std::make_pair(*vertex_it, *vertex_it2)].push_back(it->getFacetId());
		}
		else
		{
			edgesMap.insert(std::pair<Generic_Edge, std::vector<UInt> >(std::make_pair(*vertex_it2, *vertex_it), std::vector<UInt>()) );
			edgesMap[std::make_pair(*vertex_it2, *vertex_it)].push_back(it->getFacetId());
		}
	}

	for(std::map<Generic_Edge, std::vector<UInt> >::const_iterator it = edgesMap.begin(); it != edgesMap.end(); ++it)
	{
		M_edges.push_back(Edge(it->first, this, edgeId));
		M_borderEdges.push_back(Border_Edge(edgeId, this));
		M_edges.rbegin()->M_separatedFacetsIds = it->second;
		edgeId++;
	}

	edgesMap.clear();
	std::map<Generic_Edge, std::vector< std::pair<UInt, bool> > > edgesFMap;

	for(auto it = M_fractureFacets.begin(); it != M_fractureFacets.end(); ++it)
	{
		auto vertex_it  = it->getFacet().getVertexesIds().begin();
		auto vertex_it2 = vertex_it;
		++vertex_it2;

		for( ; vertex_it2 != it->getFacet().getVertexesIds().end(); ++vertex_it, ++vertex_it2)
		{
			if (*(vertex_it) < *(vertex_it2))
			{
				edgesFMap.insert(std::pair<Generic_Edge, std::vector< std::pair<UInt, bool> > >
								(std::make_pair(*vertex_it, *vertex_it2), std::vector<std::pair<UInt, bool> >())
								);
				edgesFMap[std::make_pair(*vertex_it, *vertex_it2)].push_back(std::make_pair(it->getFacetId(), it->));
			}
			else
			{
				edgesFMap.insert(std::pair<Generic_Edge, std::vector<UInt> >(std::make_pair(*vertex_it2, *vertex_it), std::vector<UInt>()) );
				edgesFMap[std::make_pair(*vertex_it2, *vertex_it)].push_back(it->getFacetId());
			}
		}

		vertex_it2 = it->getFacet().getVertexesIds().begin();

		if (*(vertex_it) < *(vertex_it2))
		{
			edgesFMap.insert(std::pair<Generic_Edge, std::vector<UInt> >(std::make_pair(*vertex_it, *vertex_it2), std::vector<UInt>()) );
			edgesFMap[std::make_pair(*vertex_it, *vertex_it2)].push_back(it->getFacetId());
		}
		else
		{
			edgesFMap.insert(std::pair<Generic_Edge, std::vector<UInt> >(std::make_pair(*vertex_it2, *vertex_it), std::vector<UInt>()) );
			edgesFMap[std::make_pair(*vertex_it2, *vertex_it)].push_back(it->getFacetId());
		}
	}

	for(std::map<Generic_Edge, std::vector<UInt> >::const_iterator it = edgesMap.begin(); it != edgesMap.end(); ++it)
	{
		M_edges.push_back(Edge(it->first, this, edgeId));
		for(std::set<UInt>::const_iterator sepIt = it->second.begin(); sepIt != it->second.end(); ++sepIt)
			if(M_facets[*sepIt].)
			sepFacets.
		M_fractureEdges.push_back(Fracture_Edge(edgeId, this));
		M_edges.rbegin()->M_separatedFacetsIds = it->second;
		edgeId++;
		sepId.clear();
	}
	*/
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

void Rigid_Mesh::showPoint(const Point & generic_point, std::ostream & out) const
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
	out << "Number of Border-Edges: " << M_borderEdges.size() <<std::endl;
	out << "Number of Fracture-Edges: " << M_fractureEdges.size() <<std::endl;
	out << "Number of Border-Fracture-Edges: " << M_borderFractureEdges.size() <<std::endl;

}

const std::vector<Rigid_Mesh::Generic_Point> Rigid_Mesh::IdToPoints(const std::vector<UInt> & pointsIds)
{
	std::vector<Generic_Point> points;
	for (auto iter : pointsIds)
		points.emplace_back(M_nodes[iter]);
	return points;
}

// --------------------   Class Edge   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Edge::Edge (const Generic_Edge & generic_edge, Geometry::Rigid_Mesh * const mesh, const UInt id):
	M_edge(generic_edge), M_mesh(mesh), M_id(id){}

Rigid_Mesh::Edge::Edge (const Edge & edge, Geometry::Rigid_Mesh * const mesh):
	M_edge(edge.getEdge()), M_mesh(mesh), M_id(edge.getId()),
	M_separatedFacetsIds(edge.getSeparatedFacetsIds()){}

Rigid_Mesh::Edge::Edge (const Edge & edge):
	M_edge(edge.getEdge()), M_mesh(edge.getMesh()), M_id(edge.getId()),
	M_separatedFacetsIds(edge.getSeparatedFacetsIds()){}

const Rigid_Mesh::Generic_Point Rigid_Mesh::Edge::getCentroid () const
{
	Generic_Point first (M_mesh->getNodesVector()[M_edge.first]);
	Generic_Point second (M_mesh->getNodesVector()[M_edge.second]);
	return (first+second) / 2.;
}

Real Rigid_Mesh::Edge::length() const
{
	Generic_Point first (M_mesh->getNodesVector()[M_edge.first]);
	Generic_Point second (M_mesh->getNodesVector()[M_edge.second]);
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
	out << "M_separatedFacets : size = " << M_separatedFacetsIds.size();
	out << "    [ ";
	if( !M_separatedFacetsIds.empty() )
		for( auto it = M_separatedFacetsIds.begin(); it != M_separatedFacetsIds.end(); ++it )
			out << *it << " ";
	out << "] " << std::endl;
}

// --------------------   Class Facet   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Facet::Facet(const Generic_Facet & generic_facet, Geometry::Rigid_Mesh * const mesh, const std::map<UInt,UInt> & old_to_new_map, const UInt m_id):
	M_mesh(mesh), M_id(m_id), M_vertexesIds(generic_facet.getVertexesVector()), M_area(generic_facet.area()),
	M_centroid(generic_facet.getCentroid()), M_unsignedNormal(generic_facet.computeNormal()), M_isFracture(generic_facet.isFracture())
{
	if(!M_mesh->M_renumber)
		for(auto it = generic_facet.getSeparatedCells().begin(); it != generic_facet.getSeparatedCells().end(); ++it)
			M_separatedCellsIds.emplace_back(*it);
	else
		for(auto it = generic_facet.getSeparatedCells().begin(); it != generic_facet.getSeparatedCells().end(); ++it)
			M_separatedCellsIds.emplace_back(old_to_new_map.at(*it));
}

Rigid_Mesh::Facet::Facet(const Facet & facet):
	M_mesh(facet.getMesh()), M_id(facet.getId()), M_vertexesIds(facet.getVertexesIds()),
	M_separatedCellsIds(facet.getSeparatedCellsIds()), M_area(facet.area()), M_centroid(facet.getCentroid()),
	M_unsignedNormal(facet.getUnsignedNormal()), M_isFracture(facet.isFracture()){}

Rigid_Mesh::Facet::Facet(const Facet & facet, Geometry::Rigid_Mesh * const mesh):
	M_mesh(mesh), M_id(facet.getId()), M_vertexesIds(facet.getVertexesIds()),
	M_separatedCellsIds(facet.getSeparatedCellsIds()), M_area(facet.area()), M_centroid(facet.getCentroid()),
	M_unsignedNormal(facet.getUnsignedNormal()), M_isFracture(facet.isFracture()){}

// ==================================================
// Methods
// ==================================================

void Rigid_Mesh::Facet::showMe(std::ostream & out) const
{
	out << "Type = Facet: " << std::endl;
	out << "ID: " << M_id << std::endl;
	out << "M_mesh = " << M_mesh << std::endl;
	out << "Number of Nodes: " << M_vertexesIds.size()<<" : "<<std::endl;
	UInt counter=1;
	for(auto j = M_vertexesIds.begin(); j!=M_vertexesIds.end(); ++j)
	{
		out << "M_idP"<<counter<<" = " <<*(j)<<"\t";
		++counter;
	}
	out << std::endl;
	out << "Size: " << M_area << std::endl;
	out << "Center: ";
	M_mesh->showPoint(M_centroid, out);
	out << "M_separatedCells : size = " << M_separatedCellsIds.size();
	out << "    [ ";
	if( !M_separatedCellsIds.empty() )
		for( auto it = M_separatedCellsIds.begin(); it != M_separatedCellsIds.end(); ++it )
			out << *it << " ";
	out << "] " << std::endl;
}

// --------------------   Class Cell   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Cell::Cell( const Generic_Cell & generic_cell, Geometry::Rigid_Mesh * const mesh, const UInt m_id):
	M_mesh(mesh), M_Id(m_id), M_zoneCode(generic_cell.getZoneCode()),
	M_vertexesIds(generic_cell.getVertexesVector()),
	M_centroid(generic_cell.getCentroid()),
	M_volume(generic_cell.volume())
{
	for(auto it = generic_cell.getFacetsSet().begin(); it != generic_cell.getFacetsSet().end(); ++it)
		M_facetsIds.emplace_back(*it);
	for(auto it = generic_cell.getNeighborsSet().begin(); it != generic_cell.getNeighborsSet().end(); ++it)
		M_neighborsIds.emplace_back(*it);
}

Rigid_Mesh::Cell::Cell(const Cell & cell, Geometry::Rigid_Mesh * const mesh):
	M_mesh(mesh), M_Id(cell.getId()), M_zoneCode(cell.getZoneCode()),
	M_vertexesIds(cell.getVertexesIds()), M_facetsIds(cell.getFacetsIds()),
	M_neighborsIds(cell.getNeighborsIds()), M_centroid(cell.getCentroid()),
	M_volume(cell.getVolume()) {}

Rigid_Mesh::Cell::Cell(const Cell & cell):
	M_mesh(cell.getMesh()), M_Id(cell.getId()), M_zoneCode(cell.getZoneCode()),
	M_vertexesIds(cell.getVertexesIds()), M_facetsIds(cell.getFacetsIds()),
	M_neighborsIds(cell.getNeighborsIds()), M_centroid(cell.getCentroid()),
	M_volume(cell.getVolume()) {}

// ==================================================
// Methods
// ==================================================

void Rigid_Mesh::Cell::showMe(std::ostream  & out) const
{
	out << "Type = Cell :" << std::endl;
	out << "ID: " << M_Id << std::endl;
	out << "M_mesh = " << M_mesh << std::endl;
	out << "Number of Nodes: " << M_vertexesIds.size()<<" : "<<std::endl;
	UInt counter=1;
	for(auto j = M_vertexesIds.begin(); j!=M_vertexesIds.end(); ++j)
	{
		out << "M_idP"<<counter<<" = " <<*(j)<<"\t";
		++counter;
	}
	out << std::endl;
	out << "Volume: " << M_volume << std::endl;
	out << "Centroid: ";
	M_mesh->showPoint(M_centroid, out);
	out << "M_Neighbors_Ids : size = " << M_neighborsIds.size();
	out << "    [ ";
	if( !M_neighborsIds.empty() )
		for( auto it = M_neighborsIds.begin(); it != M_neighborsIds.end(); ++it )
			out << *it << " ";
	out << "] " << std::endl;
}

// --------------------   Class Edge_ID   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Edge_ID::Edge_ID(const UInt edge_id, Geometry::Rigid_Mesh * const mesh):
	Edge_Id(edge_id), M_mesh(mesh){}

// --------------------   Class Border_Edge   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Border_Edge::Border_Edge (const UInt edge_Id, Geometry::Rigid_Mesh * const mesh):
	Edge_ID(edge_Id, mesh){}

Rigid_Mesh::Border_Edge::Border_Edge (const Border_Edge & border_edge, Geometry::Rigid_Mesh * const mesh):
	Edge_ID(border_edge.getEdgeId(), mesh){}

Rigid_Mesh::Border_Edge::Border_Edge (const Border_Edge & e):
	Edge_ID(e.getEdgeId(), e.getMesh()){}

// --------------------   Class Fracture_Edge   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Fracture_Edge::Fracture_Edge (const UInt edge_Id, const std::set<UInt> & fracture_Ids, Geometry::Rigid_Mesh * const mesh):
	Edge_ID(edge_Id, mesh)
{
	for (auto it = fracture_Ids.begin(); it != fracture_Ids.end(); ++it)
		Fracture_Ids.emplace_back(*(it));
}

Rigid_Mesh::Fracture_Edge::Fracture_Edge (const Fracture_Edge & fracture_edge, Geometry::Rigid_Mesh * const mesh):
	Edge_ID(fracture_edge.getEdgeId(), mesh), Fracture_Ids(fracture_edge.getFractureIds()){}

Rigid_Mesh::Fracture_Edge::Fracture_Edge (const Fracture_Edge & e):
	Edge_ID(e.getEdgeId(), e.getMesh()), Fracture_Ids(e.getFractureIds()){}

// --------------------   Class Border_Fracture_Edge   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Border_Fracture_Edge::Border_Fracture_Edge (const UInt edge_Id, Geometry::Rigid_Mesh * const mesh):
	Edge_ID(edge_Id, mesh), M_borderEdge(nullptr), M_fractureEdge(nullptr){}

Rigid_Mesh::Border_Fracture_Edge::Border_Fracture_Edge (const UInt edge_Id, const Border_Edge & borderEdge, const Fracture_Edge & fractureEdge, Geometry::Rigid_Mesh * const mesh):
	Edge_ID(edge_Id, mesh), M_borderEdge(&borderEdge), M_fractureEdge(&fractureEdge){}

Rigid_Mesh::Border_Fracture_Edge::Border_Fracture_Edge (const Border_Fracture_Edge & border_fracture_edge, Geometry::Rigid_Mesh * const mesh):
	Edge_ID(border_fracture_edge.getEdgeId(), mesh),
	M_borderEdge(&(border_fracture_edge.getBorderEdge())),
	M_fractureEdge(&(border_fracture_edge.getFractureEdge())){}

Rigid_Mesh::Border_Fracture_Edge::Border_Fracture_Edge (const Border_Fracture_Edge & e):
	Edge_ID(e.getEdgeId(), e.getMesh()),
	M_borderEdge(&(e.getBorderEdge())),
	M_fractureEdge(&(e.getFractureEdge())){}

// --------------------   Class Facet_ID   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Facet_ID::Facet_ID(const UInt facet_Id, const UInt zoneCode, Geometry::Rigid_Mesh * const mesh):
	Facet_Id(facet_Id), M_zoneCode(zoneCode), M_mesh(mesh){}

// --------------------   Class Regular_Facet   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Regular_Facet::Regular_Facet(const UInt facet_Id, const UInt zoneCode, Geometry::Rigid_Mesh * const mesh):
	Facet_ID(facet_Id, zoneCode, mesh){}

Rigid_Mesh::Regular_Facet::Regular_Facet(const Regular_Facet & regular_facet):
	Facet_ID(regular_facet.getFacetId(), regular_facet.getZoneCode(), regular_facet.getMesh()) {}

Rigid_Mesh::Regular_Facet::Regular_Facet(const Regular_Facet & regular_facet, Geometry::Rigid_Mesh * const mesh):
	Facet_ID(regular_facet.getFacetId(), regular_facet.getZoneCode(), mesh) {}

// --------------------   Class Border_Facet   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Border_Facet::Border_Facet(const UInt facet_Id, const UInt zoneCode, const UInt border_Id, Geometry::Rigid_Mesh * const mesh):
	Facet_ID(facet_Id, zoneCode, mesh), Border_Id(border_Id){}

Rigid_Mesh::Border_Facet::Border_Facet(const Border_Facet & border_facet, Geometry::Rigid_Mesh * const mesh):
	Facet_ID(border_facet.getFacetId(), border_facet.getZoneCode(), mesh), Border_Id(border_facet.getBorderId()){}

Rigid_Mesh::Border_Facet::Border_Facet(const Border_Facet & border_facet):
	Facet_ID(border_facet.getFacetId(), border_facet.getZoneCode(), border_facet.getMesh()), Border_Id(border_facet.getBorderId()) {}

// --------------------   Class Fracture_Facet   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Fracture_Facet::Fracture_Facet(const std::tuple<UInt,UInt,UInt,UInt> facet_Ids, const std::set<UInt> & fracture_Ids, Geometry::Rigid_Mesh * const mesh):
	Facet_ID(std::get<0>(facet_Ids), std::get<3>(facet_Ids), mesh), Cells_number(std::get<2>(facet_Ids)), M_Id(std::get<1>(facet_Ids))
{
	for (auto it = fracture_Ids.begin(); it != fracture_Ids.end(); ++it)
		Fracture_Ids.emplace_back(*(it));
}

Rigid_Mesh::Fracture_Facet::Fracture_Facet(const Fracture_Facet & fracture_facet, Geometry::Rigid_Mesh * const mesh):
	Facet_ID(fracture_facet.getFacetId(), fracture_facet.getZoneCode(), mesh),
	Cells_number(fracture_facet.getCellsnumber()), M_Id(fracture_facet.getId()),
	Fracture_Ids(fracture_facet.getFractureIds()),
	Fracture_Neighbors(fracture_facet.getFractureNeighbors()) {}

Rigid_Mesh::Fracture_Facet::Fracture_Facet(const Fracture_Facet & fracture_facet):
	Facet_ID(fracture_facet.getFacetId(), fracture_facet.getZoneCode(), fracture_facet.getMesh()),
	Cells_number(fracture_facet.getCellsnumber()), M_Id(fracture_facet.getId()),
	Fracture_Ids(fracture_facet.getFractureIds()),
	Fracture_Neighbors(fracture_facet.getFractureNeighbors()) {}

}// namespace Geometry
