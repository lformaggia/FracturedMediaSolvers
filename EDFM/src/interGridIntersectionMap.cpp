 /*!
 *	@file interGridIntersectionMap.cpp
 *	@brief Class for storage of grid intersections without repetitions (definitions).
 *
 *	@author Luca Turconi <lturconi@gmail.com>
 *  @date 27-09-2012
 *
 */

#include <iomanip>
 
#include "interGridIntersectionMap.hpp"
 
namespace Intersect
{
	// --------------------   Class GridIntersectionMap   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================
	GridIntersectionMap::GridIntersectionMap() : M_IntersectionMap() {}
	
	GridIntersectionMap::GridIntersectionMap(const UInt & Nx, const UInt & Ny, const UInt & Nz) :
		M_Nx(Nx), M_Ny(Ny), M_Nz(Nz), M_IntersectionMap() {}
	
	GridIntersectionMap::GridIntersectionMap(const Geometry::CPgrid & grid) :
		M_Nx(grid.Nx()), M_Ny(grid.Ny()), M_Nz(grid.Nz()) {}
	
	GridIntersectionMap::GridIntersectionMap(const GridIntersectionMap & gim) :
		M_IntersectionMap(gim.getIntersectionMap()) {}
	
	GridIntersectionMap::~GridIntersectionMap() {}
	
// ==================================================
// Methods
// ==================================================
	void GridIntersectionMap::insert(const GridIntersectionMapElement & mapElem)
	{
		std::pair<GridIntersectionMap_Iterator_Type,bool> ret;
		
		ret = M_IntersectionMap.insert( mapElem );
		if (ret.second==false)
		{
			ret.first->second.merge(mapElem.second);
		}
	}
	
	bool GridIntersectionMap::exportToGridIntersections(GridIntersections & gridInter) const
	{
		if( gridInter.Nx()!=M_Nx || gridInter.Ny()!=M_Ny || gridInter.Nz()!=M_Nz )
		{
			std::cerr << " *** Error: dimensions of gridInter and GridIntersectionMap must agree! ***" << std::endl;
			return 0;
		}
		
		gridInter.clearAll();
		
		GridIntersections_Iterator_Type git;
		CellIntersections cellInter;
		Intersection inter;
		
		bool insDone;
		
		for(GridIntersectionMap_Const_Iterator_Type it = this->begin();
			it != this->end(); ++it)
		{
			for(Linker_Const_Iterator_Type lit = it->second.begin();
				lit != it->second.end(); ++lit)
			{
				// lit->first = The cell ID
				cellInter.i() = gridInter.i(lit->first);
				cellInter.j() = gridInter.j(lit->first);
				cellInter.k() = gridInter.k(lit->first);
				
				inter.first = lit->second;	// The Edge ID
				inter.second = it->first;	// The Point3D
				
				cellInter.insert(inter);
				
				insDone = gridInter.insert(cellInter,0);
				
				if(!insDone)
				{
					git = gridInter.find(lit->first);
					
					git->second.importIntersections(cellInter);
					gridInter.incrementNintersections();
				}
				
				cellInter.clearIntersections();
			}
		}
		
		return 1;
	}
	
	bool GridIntersectionMap::importFromGridIntersections(GridIntersections & gridInter)
	{
		if( gridInter.Nx()!=M_Nx || gridInter.Ny()!=M_Ny || gridInter.Nz()!=M_Nz )
		{
			std::cerr << " *** Error: dimensions of gridInter and GridIntersectionMap must agree! ***" << std::endl;
			return 0;
		}
		
		GridIntersectionMapElement elem;
		LinkerElement link;
		
		for(GridIntersections_Const_Iterator_Type git=gridInter.begin();
			git!=gridInter.end(); ++git)
		{
			for(CellIntersections_Const_Iterator_Type cit=git->second.begin();
				cit!=git->second.end(); ++cit)
			{
				//git->first = IDcell
				// cit->first = IDedge
				// cit->second = point3D
				link.first = git->first;
				link.second = cit->first;
				
				elem.first = cit->second;
				elem.second.insert(link);
				
				this->insert(elem);
				
				elem.second.clearLinksList();
			}
		}
		
		return 1;
	}
	
	bool GridIntersectionMap::exportVtk(const std::string & fileName) const
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
		
		UInt nPoints = M_IntersectionMap.size();
		UInt CellType = 1; // for VTK_POINT
		
		// Header
		filestr << "# vtk DataFile Version 3.1" << std::endl;
		filestr << "this is a file created for Paraview" << std::endl;
		filestr << "ASCII" << std::endl;
		filestr << "DATASET UNSTRUCTURED_GRID" << std::endl;
		filestr << std::endl;	// The fifth line is empty.
		
		filestr << std::scientific << std::setprecision(10);
		
		// Pointdata
		filestr << "POINTS " << nPoints << " double" << std::endl;
		
		for(GridIntersectionMap_Const_Iterator_Type it=this->begin();
			it!=this->end(); ++it)
		{
			filestr << it->first.x << " "
					<< it->first.y << " "
					<< it->first.z << std::endl;
		
		}
		
		filestr << std::endl;
		
		// Celldata
		filestr << "CELLS " << nPoints << " " << 2*nPoints << std::endl;
		
		for(UInt i=0; i<nPoints; ++i)
		{
			filestr << "1 " << i << std::endl;
		}
		
		filestr << std::endl;
		
		filestr << "CELL_TYPES " << nPoints << std::endl;
			
		for(UInt i=0; i<nPoints; ++i)
			filestr << CellType<< std::endl;
		filestr << std::endl;
		
		filestr.close();
		
		return 1;
	}	
	
	void GridIntersectionMap::showMe(std::ostream  & out) const
	{
		out << " GridIntersectionMap :" << std::endl;
		out << " grid dimensions : Nx = " << M_Nx << ", Ny = " << M_Ny;
		out << ", Nz = " << M_Nz << std::endl;
		out << " #intersections = " << this->size() << std::endl;
		
		out << std::endl;
		out << " -----------------------------------------------------" << std::endl;
		
		for(GridIntersectionMap_Const_Iterator_Type it = this->begin();
			it != this->end(); ++it)
		{
			out << " Intersection : " << it->first << std::endl;
			it->second.showMe();
		}
		
		out << " -----------------------------------------------------" << std::endl;
	}
	
} // namespace Intersect