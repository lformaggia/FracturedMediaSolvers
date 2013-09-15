 /*!
 *	@file GeomExportVtkTool.cpp
 *	@brief Some function to Export in Vtk the elementary geometrical objects.
 *
 *	@author Luca Turconi <luca.turconi@moxoff.com>
 *
 */ 

#include "GeomExportVtkTool.hpp"

#include<fstream>

namespace Geometry{

	bool exportSetOfPointsVtk( const std::vector<Geometry::Point2D> & set,
							   const std::string & fileName )
	{
		std::fstream filestr;
		
		filestr.open (fileName.c_str(), std::ios_base::out);
	
		if (filestr.is_open())
		{
			std::cout << std::endl << " File: " << fileName << ", successfully opened";
		}
		else
		{
			std::cerr << std::endl << " *** Error: file not opened *** "
					  << std::endl << std::endl;
			return  0;
		}
		
		std::cout << std::endl <<
			" Exporting set of points in Vtk format... "
			<< std::endl;
		
		UInt nPoints( set.size() );
		
		UInt nCells(nPoints);
		
		UInt CellType = 1; // for VTK_VERTEX
		
		// Header
		filestr << "# vtk DataFile Version 3.1" << std::endl;
		filestr << "this is a file created for Paraview" << std::endl;
		filestr << "ASCII" << std::endl;
		filestr << "DATASET UNSTRUCTURED_GRID" << std::endl;
		filestr << std::endl;	// The fifth line is empty.
		
		filestr << std::scientific << std::setprecision(10);
		
		// Pointdata
		filestr << "POINTS " << nPoints << " double" << std::endl;
		for( std::vector<Geometry::Point2D>::const_iterator it = set.begin();
		     it != set.end(); ++it )
		{
			filestr << it->x() << " " << it->y() << " 0" << std::endl;
			filestr << std::endl;
		}
		
		// Celldata
		UInt nVal(2*nCells);
		
		filestr << "CELLS " << nCells << " " << nVal << std::endl;
		
		UInt id(0);
		
		for( std::vector<Geometry::Point2D>::const_iterator it = set.begin();
		     it != set.end(); ++it )
		{
			filestr << "1 " << id << std::endl;
			id++;
		}
		
		filestr << std::endl;
		
		filestr << "CELL_TYPES " << nCells << std::endl;
		for( std::vector<Geometry::Point2D>::const_iterator it = set.begin();
		     it != set.end(); ++it )
		{
			filestr << CellType << " ";
		}
		filestr << std::endl;
		
		filestr.close();
		
		return 1;
	
	}	

	bool exportSegment2DVtk( const Geometry::Segment2D & seg,
							 const std::string & fileName )
	{
		std::fstream filestr;
		
		filestr.open (fileName.c_str(), std::ios_base::out);
	
		if (filestr.is_open())
		{
			std::cout << std::endl << " File: " << fileName << ", successfully opened";
		}
		else
		{
			std::cerr << std::endl << " *** Error: file not opened *** "
					  << std::endl << std::endl;
			return  0;
		}
		
		std::cout << std::endl <<
			" Exporting Geometry::Segment2D in Vtk format... "
			<< std::endl;
		
		UInt nPoints(2);
		
		UInt nCells(1);
		
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
		filestr << seg.source().x() << " " << seg.source().y()
				<< " 0" << std::endl;
		filestr << seg.target().x() << " " << seg.target().y()
				<< " 0" << std::endl;
		filestr << std::endl;
		
		// Celldata
		filestr << "CELLS " << nCells << " " << 3 << std::endl;
		filestr << "2 0 1" << std::endl;
		filestr << std::endl;
		
		filestr << "CELL_TYPES " << nCells << std::endl;
		filestr << CellType << std::endl;
		filestr << std::endl;
		
		filestr.close();
		
		return 1;
	}
	
} // namespace Geometry
