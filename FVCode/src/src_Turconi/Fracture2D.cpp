 /*!
 *	@file Fracture2D.cpp
 *	@brief Class for fracture and fracture network in 2D space (definition).
 *
 *	@author Luca Turconi <luca.turconi@moxoff.com>
 *
 */ 

#include "Fracture2D.hpp"

#include<fstream>
//#include<limits>
//#include<cstring>


namespace Geometry{

	// --------------------   Class Fracture2D   --------------------
	
// ==================================================
// Constructors & Destructor
// ==================================================
	Fracture2D::Fracture2D() {}
	
	Fracture2D::Fracture2D(const Geometry::Fracture2D & f) :
		M_segment(f.segment()), M_h(f.h()),
		M_Nelements(f.Nelements()), M_nodes(f.nodes()),
		M_permeability(f.permeability()), M_compressibility(f.compressibility()),
		M_aperture(f.aperture()) {}
	
	Fracture2D::Fracture2D(const Geometry::Point2D source, const Geometry::Point2D target) :
		M_segment(source,target) {}
	
	Fracture2D::~Fracture2D() {}
		
// ==================================================
// Methods
// ==================================================
	Real Fracture2D::buildFractureDiscretization(const Real & h)
	{
		return this->buildFractureDiscretization(
			static_cast<UInt>((std::sqrt(M_segment.squared_length())/h)+1) );
	}
	
	Real Fracture2D::buildFractureDiscretization(const UInt & Nelements)
	{
		M_Nelements = Nelements;
	
		Real length = std::sqrt(M_segment.squared_length());
		
		M_h = length / static_cast<Real>(M_Nelements);
		
		// Fill the nodes vector
		M_nodes.resize(M_Nelements+1);
		M_nodes[0] = M_segment.source();
		
//		std::cout << " length = " << length << std::endl;
//		std::cout << " dx = " << M_segment.direction().dx()/length << std::endl;
//		std::cout << " dy = " << M_segment.direction().dy()/length << std::endl;

//		std::cout << " nodes[0] = " << M_nodes[0] << std::endl;

		for(UInt i=1; i<=M_Nelements; ++i)
		{
			M_nodes[i] = Geometry::Point2D(
							M_segment.source().x() + i*M_h*M_segment.direction().dx()/length,
							M_segment.source().y() + i*M_h*M_segment.direction().dy()/length );
											
//			std::cout << " nodes[" << i << "] = " << M_nodes[i] << std::endl;

		}
		
		return M_h;
	}
	
	bool Fracture2D::exportVtk(const std::string & fileName) const
	{
		std::fstream filestr;
		
		filestr.open (fileName.c_str(), std::ios_base::out);
	
		if (filestr.is_open())
		{
			std::cout << std::endl << " File: " << fileName << ", successfully opened";
		}
		else
		{
			std::cerr << std::endl
					  << " *** Error: file not opened *** " << std::endl << std::endl;
			return  0;
		}
		
		std::cout << std::endl << " Exporting Fracture2D in Vtk format... " << std::endl;
		
		UInt nCells = 1;
		UInt nPoints = 2;
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
		filestr << M_segment.source().x() << " " << M_segment.source().y()
				<< " 0" << std::endl;
		filestr << M_segment.target().x() << " " << M_segment.target().y()
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
	
	bool Fracture2D::exportFractureDiscretizationVtk(const std::string & fileName) const
	{
		std::fstream filestr;
		
		filestr.open (fileName.c_str(), std::ios_base::out);
	
		if (filestr.is_open())
		{
			std::cout << std::endl << " File: " << fileName << ", successfully opened";
		}
		else
		{
			std::cerr << std::endl
					  << " *** Error: file not opened *** " << std::endl << std::endl;
			return  0;
		}
		
		std::cout << std::endl
				  << " Exporting Fracture2D discretization in Vtk format... " << std::endl;
		
		UInt nPoints = M_Nelements+1;
		UInt nCells = nPoints;
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
		for(std::vector<Geometry::Point2D>::const_iterator it = M_nodes.begin();
			it != M_nodes.end(); ++it)
		{
			filestr << it->x() << " " << it->y() << " 0" << std::endl;
		}
		filestr << std::endl;
		
		// Celldata
		filestr << "CELLS " << nCells << " " << 2*nCells << std::endl;
		for(UInt i=0; i<nPoints; ++i)
			filestr << "1 " << i << std::endl;

		filestr << std::endl;
		
		filestr << "CELL_TYPES " << nCells << std::endl;
		
		for(UInt i=1; i<=nCells; ++i)
			filestr << CellType << std::endl;
		
		filestr << std::endl;
		
		filestr.close();
		
		return 1;
	}
	
	void Fracture2D::showMe(std::ostream  & out) const
	{
		out << "Type = Fracture2D " << std::endl;
		out << " M_segment :" << std::endl;
		out << " [ ";
		out << "(" << M_segment.source().x() << "," << M_segment.source().y() << ")";
		out << " ; ";
		out << "(" << M_segment.target().x() << "," << M_segment.target().y() << ")";
		out << std::endl;
		out << " Physical Properties : " << std::endl;
		out << " -> permeability : " << M_permeability << std::endl;
		out << " -> compressibility : " << M_compressibility << std::endl;
		out << " -> aperture : " << M_aperture << std::endl;
	}
	
} // namespace Geometry
