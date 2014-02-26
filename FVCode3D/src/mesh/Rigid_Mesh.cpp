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
			it->second.computeCentroid();
			M_cells.push_back(Cell(it->second, this, it->first));
		}
	}
	else
	{
		UInt position = 0;
		for (auto it = generic_mesh.getCellsMap().begin(); it != generic_mesh.getCellsMap().end(); ++it)
		{
			it->second.computeCentroid();
			old_to_new_map[it->first] = position;
			M_cells.push_back(Cell(it->second, this, position));
			++position;
		}
	}

}

void Rigid_Mesh::AdjustCellNeighbors ( const std::map<UInt, UInt> & old_to_new_map )
{
	for (auto it = M_cells.begin(); it != M_cells.end(); ++it)
		for (auto Neighbors_it = it->M_Neighbors_Ids.begin(); Neighbors_it != it->M_Neighbors_Ids.end(); ++ Neighbors_it)
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

		++vertex_it;
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
				auto find_it = std::find(it->Fracture_Neighbors[m_junct].begin(),
						it->Fracture_Neighbors[m_junct].end(), it->getId());
				it->Fracture_Neighbors[m_junct].erase(find_it);
			}
		}

		++vertex_it;
		vertex_it2 = it->getFacet().getVertexesIds().begin();

		//The if state is necessary in order to order the vertexes!!
		if (*(vertex_it) < *(vertex_it2))
			m_junct = std::make_pair(*(vertex_it), *(vertex_it2));
		else
			m_junct = std::make_pair(*(vertex_it2), *(vertex_it));
		if (nodes_fracture_map.at(m_junct).size() > 1)
		{
			it->Fracture_Neighbors[m_junct] = nodes_fracture_map.at(m_junct);
			auto find_it = std::find(it->Fracture_Neighbors[m_junct].begin(),
					it->Fracture_Neighbors[m_junct].end(), it->getId());
			it->Fracture_Neighbors[m_junct].erase(find_it);
		}
	}
}

void Rigid_Mesh::M_facetsVectorsBuilder(const std::map<UInt,UInt> & old_to_new_mapCells, std::map<UInt,UInt> & old_to_new_mapFacets, Generic_Mesh & generic_mesh)
{
	UInt fractureId = 0;

	const Generic_FacetsContainer & facet_container = generic_mesh.getFacetsMap();

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
					Fracture_Facet(fracture_Ids, it->second.getRepresentedFractureIds(), this, generic_mesh)
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
					Fracture_Facet( fracture_Ids, it->second.getRepresentedFractureIds(), this, generic_mesh)
				);
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
		for (auto Facets_it = it->M_Facets_Ids.begin(); Facets_it != it->M_Facets_Ids.end(); ++ Facets_it)
			*(Facets_it) = old_to_new_mapFacets.at(*(Facets_it));
}

bool Rigid_Mesh::hasNeighborsThroughFacet(const UInt & facet_Id, const UInt & idNeighbor) const
{
	if( M_facets[facet_Id].getSeparatedCellsIds()[0] == idNeighbor ||
		M_facets[facet_Id].getSeparatedCellsIds()[1] == idNeighbor)
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
	out << "Type = Mesh3D :"<< std::endl;
	out << "code = " << this << std::endl;
	out << "Number of Nodes: " << M_nodes.size() <<std::endl;
	out << "Number of Cells: " << M_cells.size() <<std::endl;
	out << "Number of Facets: " << M_facets.size() <<std::endl;
	out << "Number of Border-Facets: " << M_borderFacets.size() <<std::endl;
	out << "Number of Fracture-Facets: " << M_fractureFacets.size() <<std::endl;
	out << "Number of Standard-Facets: " << M_internalFacets.size() <<std::endl;

}

const std::vector<Rigid_Mesh::Generic_Point> Rigid_Mesh::IdToPoints(const std::vector<UInt> & pointsIds)
{
	std::vector<Generic_Point> points;
	for (auto iter : pointsIds)
		points.emplace_back(M_nodes[iter]);
	return points;
}

// --------------------   Class Facet   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Facet::Facet(const Generic_Facet & generic_facet, Geometry::Rigid_Mesh * const mesh, const std::map<UInt,UInt> & old_to_new_map, const UInt m_id):
	M_mesh(mesh), M_Id(m_id), M_Vertexes_Ids(generic_facet.getVertexesVector()), M_size(generic_facet.area())
{
	if(!M_mesh->M_renumber)
		for(auto it = generic_facet.getSeparatedCells().begin(); it != generic_facet.getSeparatedCells().end(); ++it)
			M_separatedCells_Ids.emplace_back(*it);
	else
		for(auto it = generic_facet.getSeparatedCells().begin(); it != generic_facet.getSeparatedCells().end(); ++it)
			M_separatedCells_Ids.emplace_back(old_to_new_map.at(*it));
	computeCenter();
	computeUnsignedNormal();
}

Rigid_Mesh::Facet::Facet(const Facet & facet):
	M_mesh(facet.getMesh()), M_Id(facet.getId()), M_Vertexes_Ids(facet.getVertexesIds()),
	M_separatedCells_Ids(facet.getSeparatedCellsIds()), M_size(facet.size()), M_center(facet.getCenter()),
	M_UnsignedNormal(facet.getUnsignedNormal()){}

Rigid_Mesh::Facet::Facet(const Facet & facet, Geometry::Rigid_Mesh * const mesh):
	M_mesh(mesh), M_Id(facet.getId()), M_Vertexes_Ids(facet.getVertexesIds()),
	M_separatedCells_Ids(facet.getSeparatedCellsIds()), M_size(facet.size()), M_center(facet.getCenter()),
	M_UnsignedNormal(facet.getUnsignedNormal()){}


// ==================================================
// Protected Methods
// ==================================================

void Rigid_Mesh::Facet::computeCenter()
{
	UInt N = M_Vertexes_Ids.size();
	Generic_Point sum(0. ,0. ,0.);
	for(UInt j = 0; j < N; ++j)
		sum = sum + (M_mesh->getNodesVector()[M_Vertexes_Ids[j]])/N;
	M_center =  sum;
}

void Rigid_Mesh::Facet::computeUnsignedNormal()
{
	Generic_Vector tangent_1, tangent_2;
	tangent_1 = M_mesh->getNodesVector()[M_Vertexes_Ids[1]] - M_mesh->getNodesVector()[M_Vertexes_Ids[0]];
	tangent_2 = M_mesh->getNodesVector()[M_Vertexes_Ids[2]] - M_mesh->getNodesVector()[M_Vertexes_Ids[0]];

	Generic_Vector normal( crossProduct(tangent_1, tangent_2) );
	normal.normalize();

	M_UnsignedNormal = normal;
}

// ==================================================
// Methods
// ==================================================

void Rigid_Mesh::Facet::showMe(std::ostream & out) const
{
	out << "Type = Facet3D: " << std::endl;
	out << "ID: " << M_Id << std::endl;
	out << "M_mesh = " << M_mesh << std::endl;
	out << "Number of Nodes: " << M_Vertexes_Ids.size()<<" : "<<std::endl;
	UInt counter=1;
	for(auto j = M_Vertexes_Ids.begin(); j!=M_Vertexes_Ids.end(); ++j)
	{
		out << "M_idP"<<counter<<" = " <<*(j)<<"\t";
		++counter;
	}
	out << std::endl;
	out << "Size: " << M_size << std::endl;
	out << "Center: ";
	M_mesh->showPoint(M_center, out);
	out << "M_separatedCells : size = " << M_separatedCells_Ids.size();
	out << "    [ ";
	if( !M_separatedCells_Ids.empty() )
		for( auto it = M_separatedCells_Ids.begin(); it != M_separatedCells_Ids.end(); ++it )
			out << *it << " ";
	out << "] " << std::endl;
}

// --------------------   Class Cell   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================

Rigid_Mesh::Cell::Cell( const Generic_Cell & generic_cell, Geometry::Rigid_Mesh * const mesh, const UInt m_id):
	M_mesh(mesh), M_Id(m_id), M_zoneCode(generic_cell.getZoneCode()),
	M_Vertexes_Ids(generic_cell.getVertexesVector()), M_volume(0.)
{
	for(auto it = generic_cell.getFacetsSet().begin(); it != generic_cell.getFacetsSet().end(); ++it)
		M_Facets_Ids.emplace_back(*it);
	for(auto it = generic_cell.getNeighborsSet().begin(); it != generic_cell.getNeighborsSet().end(); ++it)
		M_Neighbors_Ids.emplace_back(*it);
	M_centroid = generic_cell.getCentroid();
	M_volume = generic_cell.volume();
}

Rigid_Mesh::Cell::Cell(const Cell & cell, Geometry::Rigid_Mesh * const mesh):
	M_mesh(mesh), M_Id(cell.getId()), M_Vertexes_Ids(cell.getVertexesIds()),
	M_Neighbors_Ids(cell.getNeighborsIds()), M_Facets_Ids(cell.getFacetsIds()),
	M_centroid(cell.getCentroid()), M_volume(cell.getVolume()), M_zoneCode(cell.getZoneCode())
{}

Rigid_Mesh::Cell::Cell(const Cell & cell):
	M_mesh(cell.getMesh()), M_Id(cell.getId()), M_Vertexes_Ids(cell.getVertexesIds()),
	M_Neighbors_Ids(cell.getNeighborsIds()), M_Facets_Ids(cell.getFacetsIds()),
	M_centroid(cell.getCentroid()), M_volume(cell.getVolume()), M_zoneCode(cell.getZoneCode())
{}

// ==================================================
// Methods
// ==================================================

void Rigid_Mesh::Cell::showMe(std::ostream  & out) const
{
	out << "Type = Cell3D :" << std::endl;
	out << "ID: " << M_Id << std::endl;
	out << "M_mesh = " << M_mesh << std::endl;
	out << "Number of Nodes: " << M_Vertexes_Ids.size()<<" : "<<std::endl;
	UInt counter=1;
	for(auto j = M_Vertexes_Ids.begin(); j!=M_Vertexes_Ids.end(); ++j)
	{
		out << "M_idP"<<counter<<" = " <<*(j)<<"\t";
		++counter;
	}
	out << std::endl;
	out << "Volume: " << M_volume << std::endl;
	out << "Centroid: ";
	M_mesh->showPoint(M_centroid, out);
	out << "M_Neighbors_Ids : size = " << M_Neighbors_Ids.size();
	out << "    [ ";
	if( !M_Neighbors_Ids.empty() )
		for( auto it = M_Neighbors_Ids.begin(); it != M_Neighbors_Ids.end(); ++it )
			out << *it << " ";
	out << "] " << std::endl;
}

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

Rigid_Mesh::Fracture_Facet::Fracture_Facet(const std::tuple<UInt,UInt,UInt,UInt> facet_Ids, const std::set<UInt> & fracture_Ids, Geometry::Rigid_Mesh * const mesh, Generic_Mesh & generic_mesh):
	Facet_ID(std::get<0>(facet_Ids), std::get<3>(facet_Ids), mesh), Cells_number(std::get<2>(facet_Ids)), M_Id(std::get<1>(facet_Ids))
{
	for (auto it = fracture_Ids.begin(); it != fracture_Ids.end(); ++it)
		Fracture_Ids.emplace_back(*(it));
}

Rigid_Mesh::Fracture_Facet::Fracture_Facet(const Fracture_Facet & fracture_facet, Geometry::Rigid_Mesh * const mesh):
	Facet_ID(fracture_facet.getFacetId(), fracture_facet.getZoneCode(), mesh), Cells_number(fracture_facet.getCellsnumber()),
	M_Id(fracture_facet.getId()),
	Fracture_Ids(fracture_facet.getFractureIds()), Fracture_Neighbors(fracture_facet.getFractureNeighbors()) {}

Rigid_Mesh::Fracture_Facet::Fracture_Facet(const Fracture_Facet & fracture_facet):
	Facet_ID(fracture_facet.getFacetId(), fracture_facet.getZoneCode(), fracture_facet.getMesh()), Cells_number(fracture_facet.getCellsnumber()),
	M_Id(fracture_facet.getId()),
	Fracture_Ids(fracture_facet.getFractureIds()), Fracture_Neighbors(fracture_facet.getFractureNeighbors()) {}

}// namespace Geometry