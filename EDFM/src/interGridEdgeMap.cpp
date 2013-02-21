 /*!
 *	@file interGridEdgeMap.cpp
 *	@brief Class for storage of grid edges without repetitions (definition).
 *
 *	@author Luca Turconi <lturconi@gmail.com>
 *  @date 27-09-2012
 *
 */

#include <iomanip>
 
#include "interGridEdgeMap.hpp"

namespace Intersect
{
	// --------------------   Class Linker   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================
	Linker::Linker() : M_LinksList() {}
	
	Linker::Linker(const Linker & l) : M_LinksList(l.getLinksList()) {}
	
	Linker::~Linker() {}
	
// ==================================================
// Methods
// ==================================================
	void Linker::merge(const Linker & l)
	{
		M_LinksList.insert(this->begin(),l.begin(),l.end());
	}
	
	void Linker::showMe(std::ostream  & out) const
	{
		out << " LINKER: " << std::endl;
		for(Linker_Const_Iterator_Type it = this->begin();
			it != this->end(); ++it)
		{
			out << " idCell = " << it->first << " idEdge = " << it->second << std::endl;
		}
	}

	// --------------------   Class GridEdgeMap   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================
	GridEdgeMap::GridEdgeMap() : M_EdgeMap() {}
	
	GridEdgeMap::GridEdgeMap(const UInt & Nx, const UInt & Ny, const UInt & Nz) :
		M_Nx(Nx), M_Ny(Ny), M_Nz(Nz), M_EdgeMap() {}
	
	GridEdgeMap::GridEdgeMap(const Geometry::CPgrid & grid) :
		M_Nx(grid.Nx()), M_Ny(grid.Ny()), M_Nz(grid.Nz())
		{ this->extractEdges(grid); }
	
	GridEdgeMap::GridEdgeMap(const GridEdgeMap & gem) :
		M_EdgeMap(gem.getEdgeMap()) {}
	
	GridEdgeMap::~GridEdgeMap() {}
	
// ==================================================
// Methods
// ==================================================
	bool GridEdgeMap::extractEdges(const Geometry::CPgrid & g)
	{
		if( g.Nx()!=M_Nx || g.Ny()!=M_Ny || g.Nz()!=M_Nz )
		{
			std::cerr << " *** Error: dimensions of grid and GridEdgeMap must agree! ***" << std::endl;
			return 0;
		}
		
		Geometry::CPcell cell;
		GridEdgeMapElement GEMelem;
		LinkerElement link;
		UInt idCell;
		
		for(UInt i=1; i<=g.Nx(); ++i)
		{
			for(UInt j=1; j<=g.Ny(); ++j)
			{
				for(UInt k=1; k<=g.Nz(); ++k)
				{
					if( g.cell(i,j,k).getActnum() )
					{
						idCell = i + (j-1)*M_Nx + (k-1)*M_Nx*M_Ny;
						cell = g.cell(i,j,k);
						for(UInt idEdge=1; idEdge<=12; ++idEdge)
						{
							GEMelem.first = cell.getEdge(idEdge);
							link.first = idCell;
							link.second = idEdge;
							GEMelem.second.insert(link);
						
							this->insert(GEMelem);
						
							GEMelem.second.clearLinksList();
						}
					}
				}
			}
		}
		
		return 1;
	}
	
	void GridEdgeMap::insert(const GridEdgeMapElement & mapElem)
	{
		std::pair<GridEdgeMap_Iterator_Type,bool> ret;
		
		ret = M_EdgeMap.insert( mapElem );
		if (ret.second==false)
		{
			ret.first->second.merge(mapElem.second);
		}
	}
	
	bool GridEdgeMap::exportVtk(const std::string & fileName) const
	{
		std::fstream filestr;
		
		filestr.open (fileName.c_str(), std::ios_base::out);
	
		if (filestr.is_open())
		{
			std::cout << std::endl << " File: " << fileName << ", successfully opened";
		}
		else
		{
			std::cerr << std::endl << " *** Error: file not opened *** " << std::endl << std::endl;
			return  0;
		}
		
		std::cout << std::endl << " Exporting grid in Vtk format... " << std::endl;
		
		UInt nPoints = 2*M_EdgeMap.size();
		UInt nCells = M_EdgeMap.size();
		UInt CellType = 3; // for VTK_LINE
		
		// Header
		filestr << "# vtk DataFile Version 3.1" << std::endl;
		filestr << "this is a file created for Paraview" << std::endl;
		filestr << "ASCII" << std::endl;
		filestr << "DATASET UNSTRUCTURED_GRID" << std::endl;
		filestr << std::endl;	// The fifth line is empty.
		
		filestr << std::scientific << std::setprecision(10);
		
		// Pointdata
		filestr << "POINTS " << nPoints << " double" << std::endl;
		
		for(GridEdgeMap_Const_Iterator_Type it=this->begin();
			it!=this->end(); ++it)
		{
			filestr << it->first.A().x << " "
					<< it->first.A().y << " "
					<< it->first.A().z << std::endl;
		
			filestr << it->first.B().x << " "
					<< it->first.B().y << " "
					<< it->first.B().z << std::endl;
		}
		
		filestr << std::endl;
		
		// Celldata
		filestr << "CELLS " << nCells << " " << 3*nCells << std::endl;
		
		for(UInt i=0; i<nCells; ++i)
		{
			filestr << "2 " << 2*i << " " << 2*i+1 << std::endl;
		}
		
		filestr << std::endl;
		
		filestr << "CELL_TYPES " << nCells << std::endl;
			
		for(UInt i=0; i<nCells; ++i)
			filestr << CellType<< std::endl;
		filestr << std::endl;
		
		filestr.close();
		
		return 1;
	}
	
	void GridEdgeMap::showMe(std::ostream  & out) const
	{
		out << " GridEdgeMap :" << std::endl;
		out << " grid dimensions : Nx = " << M_Nx << ", Ny = " << M_Ny;
		out << ", Nz = " << M_Nz << std::endl;
		out << " #edges = " << this->size() << std::endl;
		
		out << std::endl;
		out << " -----------------------------------------------------" << std::endl;
		
		for(GridEdgeMap_Const_Iterator_Type it = this->begin();
			it != this->end(); ++it)
		{
			out << " Edge : " << it->first << std::endl;
			it->second.showMe();
		}
		
		out << " -----------------------------------------------------" << std::endl;
	}
	
} // namespace Intersect
