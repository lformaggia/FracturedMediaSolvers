 /*!
 *	@file Fracture2D.cpp
 *	@brief Class for fracture network in 2D space (definition).
 *
 *	@author Luca Turconi <luca.turconi@moxoff.com>
 *
 */ 
 
#include "FractureNetwork2D.hpp"

#include <boost/lexical_cast.hpp>

#include<fstream>
#include<limits>
#include<cstring>

 
namespace Geometry{
 
 // --------------------   Class FractureNetwork2D   --------------------
	
// ==================================================
// Constructors & Destructor
// ==================================================
	FractureNetwork2D::FractureNetwork2D() {}
	
	FractureNetwork2D::FractureNetwork2D(const Geometry::FractureNetwork2D & fn) :
		M_fractureNetwork(fn.getNetwork()) {}
	
	FractureNetwork2D::FractureNetwork2D(const std::string & fabFileName)
	{
		// create the input file stream
		std::ifstream myfile;
		std::string line;
		bool flag(0);
		
		// open .fab file
		myfile.open(fabFileName.c_str());
		
		// extract the first two line
		getline (myfile,line);
		getline (myfile,line);
		
		// extract the third line (the first interesting thing is the unit system used)
		getline (myfile,line);

		// unit system identification
		std::string::size_type pos11;
		pos11 = line.find_first_of("=");
	    char provv11[1];
		line.copy(provv11, 1, pos11+2);
		if( std::strcmp(provv11, "M") == 0 )
		{
			M_isMetric = 1;
			std::cout << "METRICO" << std::endl;
		}
		
		UInt nfractures;
		
		while( myfile.good() && !flag ){
	
			getline (myfile,line);
			std::string::size_type pos, pos1;
			
 			// search the "No_Fractures" attribute
			pos = line.find("No_Fractures");

 			// if find the "No_Fractures" attribute
 			if( pos != std::string::npos )
 			{
 				// extract "No_Fractures" attribute
				pos = line.find_first_of("0123456789");
				pos1 = line.find_last_of("0123456789");
				
				char provv[pos1+2-pos];
				std::string::size_type length = line.copy(provv, pos1+1-pos , pos);
				provv[length] = '\0';
				nfractures = atoi(provv);
			}
			
 			// search the "BEGIN FRACTURE" section
			pos = line.find("BEGIN FRACTURE");
			// if find the "BEGIN FRACTURE" section
 			if( pos != std::string::npos )
 			{
 				// set the flag value to 1
//				std::cout <<"finitoooo"<<std::endl;
				flag = 1;
			}
		}
		
		// when it exits from while loop
		// we are at the beginning of section "BEGIN FRACTURE"
		
		for( UInt i=0; i<nfractures; ++i )
		{
			getline (myfile,line);
			UInt n_points;
			Real perm, compr, aperture;
			FabFile::readProperties(line, n_points, perm, compr, aperture);
			std::vector<Geometry::Point2D> points;
			for( UInt j=0; j<n_points; ++j )
			{
				getline (myfile,line);
				Real px,py,pz;
				FabFile::readPoints(line,px,py,pz);
				Geometry::Point2D p(px,py);
				points.push_back(p);
//				std::cout << "p = ( " << px << " , " << py << " )" << std::endl;
			}	
			
			// get the last line of the fracture (index 0)
			// I don't know which informations it contains
			getline (myfile,line);
			
			// create the new fracture
			Fracture2D f( points[0], points[1] );
			// fill its physical properties
			f.setPermeability(perm);
			f.setCompressibility(compr);
			f.setAperture(aperture);
			
			// add the new fracture to the fracture network
			M_fractureNetwork.push_back(f);

		}
			
		// forget we hit the end of file
		myfile.clear();
		// move to the start of the file
		myfile.seekg(0, std::ios::beg);
	}
	
	FractureNetwork2D::~FractureNetwork2D() {}
	
// ==================================================
// Methods
// ==================================================
	Real FractureNetwork2D::distance(const Geometry::Point2D & p) const
	{
		Real minDistance(std::numeric_limits<Real>::max());
		
		for(std::vector<Geometry::Fracture2D>::const_iterator it = M_fractureNetwork.begin();
			it != M_fractureNetwork.end(); ++it)
		{
			if( it->distance(p) < minDistance )
				minDistance = it->distance(p);
		}
		
		return minDistance;
	}
	
	void FractureNetwork2D::buildFractureNetworkDiscretization(const Real & h)
	{
		for( std::vector<Geometry::Fracture2D>::iterator it = M_fractureNetwork.begin();
			 it != M_fractureNetwork.end(); ++it )
		{
			it->buildFractureDiscretization(h);
		}
	}
	
	
	void FractureNetwork2D::buildFractureNetworkDiscretization(const UInt & Nelements)
	{
		for( std::vector<Geometry::Fracture2D>::iterator it = M_fractureNetwork.begin();
			 it != M_fractureNetwork.end(); ++it )
		{
			it->buildFractureDiscretization(Nelements);
		}
	}
	
	
	bool FractureNetwork2D::exportVtk(const std::string & filename) const
	{
		UInt i=0;
		bool status=1;
		
		for( std::vector<Geometry::Fracture2D>::const_iterator it = M_fractureNetwork.begin();
		     it != M_fractureNetwork.end() && status != 0; ++it )
		{
			status = status && it->exportVtk( filename+"_Fracture"+
							  boost::lexical_cast<std::string>(i)+
							  ".vtk" );
			++i;
		}
		
		return status;
	}
	
	bool FractureNetwork2D::exportNetworkVtk(const std::string & fileName) const
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
				  << " Exporting FractureNetwork in Vtk format... " << std::endl;
		
		UInt nCells = this->size();
		UInt nPoints = 2 * nCells;
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
		for(std::vector<Geometry::Fracture2D>::const_iterator it = M_fractureNetwork.begin();
			it != M_fractureNetwork.end(); ++it)
		{
			filestr << it->source().x() << " " << it->source().y() << " 0" << std::endl;
			filestr << it->target().x() << " " << it->target().y() << " 0" << std::endl;
		}
		filestr << std::endl;
		
		// Celldata
		filestr << "CELLS " << nCells << " " << 3*nCells << std::endl;
		for(UInt i=0; i<nCells; ++i)
			filestr << "2 " << 2*i << " " << 2*i+1 << std::endl;

		filestr << std::endl;
		
		filestr << "CELL_TYPES " << nCells << std::endl;
		
		for(UInt i=1; i<=nCells; ++i)
			filestr << CellType << std::endl;
		
		filestr << std::endl;
		
		filestr.close();
		
		return 1;
	}
	
	bool FractureNetwork2D::exportFractureDiscretizationVtk(const std::string & filename) const
	{
		UInt i=0;
		bool status=1;
		
		for( std::vector<Geometry::Fracture2D>::const_iterator it = M_fractureNetwork.begin();
		     it != M_fractureNetwork.end() && status != 0; ++it )
		{
			status = status && it->exportFractureDiscretizationVtk
							( filename+"_DiscretizationFracture"+
							  boost::lexical_cast<std::string>(i)+
							  ".vtk" );
			++i;
		}
		
		return status;
	}
	
	void FractureNetwork2D::showMe(std::ostream  & out) const
	{
		out << "Type = FractureNetwork2D " << std::endl;
		out << std::endl;
		for( UInt i = 0; i<M_fractureNetwork.size(); ++i )
		{
			out << " f[" << i << "] : ";
			M_fractureNetwork[i].showMe(out);
			out << std::endl;
		}
		
		out << std::endl;
	}
	
} // namespace Geometry


namespace FabFile{
	
	void readProperties( std::string line, UInt & n_points,
					 Real & perm, Real & compr, Real & aperture )
	{
		// extracting aperture:

		// find the end of the last number
		std::string::size_type found = line.find_last_of("0123456789");
		// cut the finally empty part of the string
		line.resize(found+1);
		// find the beginning of the last number (aperture)
		std::string::size_type found2 = line.find_last_of(" ");
		// copy the last number in provv
		char provv[found-found2];
		std::string::size_type length = line.copy(provv, found-found2 , found2+1);
		// insert the termination character
		provv[length] = '\0';
		// convert the char number in the corresponding floating point value
		aperture = atof(provv);
		// erase the last number (aparture) from the considered line
		line.resize(found2+1);
		// -------------------------------------------

		// extracting compr:

		// now repeat the same procedure to extract the new last number (compr)
		found = line.find_last_of("0123456789");
		line.resize(found+1);
		found2 = line.find_last_of(" ");
		char provv1[found-found2];
		length = line.copy(provv1, found-found2 , found2+1);
		provv1[length] = '\0';
		// convert the char number in the corresponding floating point value
		compr = atof(provv1);
		// erase the last number (compr) from the considered line
		line.resize(found2+1);
		// -------------------------------------------

		// extracting perm:

		// now repeat the same procedure to extract the new last number (perm)
		found = line.find_last_of("0123456789");
		line.resize(found+1);
		found2 = line.find_last_of(" ");
		char provv2[found-found2];
		length = line.copy(provv2, found-found2 , found2+1);
		provv2[length] = '\0';
		// convert the char number in the corresponding floating point value
		perm = atof(provv2);
		// erase the last number (perm) from the considered line
		line.resize(found2+1);

		// find the last number
		found=line.find_last_of("0123456789");
		// cut the finally empty part of the string
		line.resize(found+1);
		// find the beginning of the last number (not interested in this value)
		found2=line.find_last_of(" ");
		// erase the last number from the considered line (not interested in this value)
		line.resize(found2+1);
		// -------------------------------------------

		// extracting n_points:

		// now repeat the same procedure to extract the new last number (n_points)
		found = line.find_last_of("0123456789");
		line.resize(found+1);
		found2 = line.find_last_of(" ");
		char provv3[found-found2];
		length = line.copy(provv3, found-found2 , found2+1);
		provv3[length] = '\0';
		// convert the char number in the corresponding integer value
		n_points = atoi(provv3);
		// erase the last number from the considered line (n_points)
		line.resize(found2+1);
		// -------------------------------------------
		
		// print the number of point for the fracture definition
//		std::cout << n_points << std::endl;
		
		// end: all interesting data are extracted
	}

	void readPoints( std::string line, Real & px, Real & py, Real & pz )
	{
		// extracting pz:

		// find the end of the last number
		std::string::size_type found = line.find_last_of("0123456789");
		// cut the finally empty part of the string
		line.resize(found+1);
		// find the beginning of the last number (pz)
		std::string::size_type found2 = line.find_last_of(" ");
		// copy the last number in provv
		char provv[found-found2];
		std::string::size_type length = line.copy(provv, found-found2 , found2+1);
		// insert the termination character
		provv[length] = '\0';
		// convert the char number in the corresponding floating point value
		pz = atof(provv);
		// erase the last number (pz) from the considered line
		line.resize(found2+1);
		// -------------------------------------------

		// extracting py:
		
		found = line.find_last_of("0123456789");
		line.resize(found+1);
		found2 = line.find_last_of(" ");
		char provv1[found-found2];
		length = line.copy(provv1, found-found2 , found2+1);
		provv1[length] = '\0';
		py = atof(provv1);
		line.resize(found2+1);
		// -------------------------------------------

		// extracting px:

		found = line.find_last_of("0123456789");
		line.resize(found+1);
		found2 = line.find_last_of(" ");
		char provv2[found-found2];
		length = line.copy(provv2, found-found2 , found2+1);
		provv2[length] = '\0';
		px = atof(provv2);
		// -------------------------------------------
	
	}

} // namespace FabFile
