/*!
 * @file export.cpp
 * @brief Classes for saving files (definitions).
 */

#include "mesh/Mesh3D.hpp"
#include "mesh/Rigid_Mesh.hpp"
#include "export/export.hpp"
#include "mesh/Properties.hpp"
#include "geometry/operations.hpp"

void ExporterVTU::exportMesh(const Mesh3D & mesh, const std::string filename)
{
	std::fstream filestr;

	filestr.open (filename.c_str(), std::ios_base::out);

	if (filestr.is_open())
	{
		std::cout << std::endl << " File: " << filename << ", successfully opened";
	}
	else
	{
		std::cerr << std::endl << " *** Error: file not opened *** " << std::endl << std::endl;
		return;
	}

	std::cout << std::endl << " Exporting Mesh3D in Vtu format... " << std::endl;

	const std::vector<Geometry::Point3D> & nodes = mesh.getNodesVector();
	const std::map<UInt, Facet3D> & facets = mesh.getFacetsMap();
	const std::map<UInt, Cell3D> & cells = mesh.getCellsMap();
	UInt nPoints = nodes.size();
	UInt nCells = cells.size();
	UInt offsets = 0, faceOffsets = 0;

	// Header
	filestr << "<?xml version=\"1.0\"?>" << std::endl;
	filestr << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl;
	filestr << "\t<UnstructuredGrid>" << std::endl;
	filestr << "\t\t<Piece NumberOfPoints=\"" << nPoints << "\" NumberOfCells=\"" << nCells << "\">" << std::endl;

	// CellData
	filestr << "\t\t\t<CellData Scalars=\"cellID\">" << std::endl;
	filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"cellID\" format=\"ascii\">" << std::endl;
	filestr << std::scientific << std::setprecision(0);
	for( std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
		filestr << it->first << std::endl;
	filestr << "\t\t\t\t</DataArray>" << std::endl;
	filestr << "\t\t\t</CellData>" << std::endl;

	filestr << std::scientific << std::setprecision(10);

	// Points
	filestr << "\t\t\t<Points>" << std::endl;
	filestr << "\t\t\t\t<DataArray type=\"" << "Float32" << "\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
	for( std::vector<Geometry::Point3D>::const_iterator it = nodes.begin(); it != nodes.end(); ++it )
		filestr << it->x() << " " << it->y() << " " << it->z() <<std::endl;
	filestr << "\t\t\t\t</DataArray>" << std::endl;
	filestr << "\t\t\t</Points>" << std::endl;

	// Cells
	filestr << "\t\t\t<Cells>" << std::endl;
	//  Connectivity
	filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
	for( std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
	{
		for( std::vector<UInt>::const_iterator jt = it->second.getVertexesVector().begin(); jt != it->second.getVertexesVector().end()-1; ++jt )
			filestr << *jt << " ";
		filestr << *(it->second.getVertexesVector().rbegin()) << std::endl;
	}
	filestr << "\t\t\t\t</DataArray>" << std::endl;

	//  Offsets
	filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << std::endl;
	for( std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
	{
		offsets += it->second.vertexesNumber();
		filestr << offsets << std::endl;
	}
	filestr << "\t\t\t\t</DataArray>" << std::endl;

	//  Types
	filestr << "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
	for( UInt i=0; i < nCells ; ++i )
		filestr << "42" << std::endl;
	filestr << "\t\t\t\t</DataArray>" << std::endl;

	//  Faces
	filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"faces\" format=\"ascii\">" << std::endl;
	for( std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
	{
		filestr << it->second.facetsNumber() << std::endl;
		for( std::set<UInt>::const_iterator jt = it->second.getFacetsSet().begin(); jt != it->second.getFacetsSet().end(); ++jt )
		{
			filestr << facets.at(*jt).getNumberOfPoints() << std::endl;
			for( std::vector<UInt>::const_iterator kt = facets.at(*jt).getVertexesVector().begin(); kt != facets.at(*jt).getVertexesVector().end()-1; ++kt )
				filestr << *kt << " ";
			filestr << *(facets.at(*jt).getVertexesVector().rbegin()) << std::endl;
		}
	}
	filestr << "\t\t\t\t</DataArray>" << std::endl;

	//  Faceoffsets
	filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"faceoffsets\" format=\"ascii\">" << std::endl;
	for( std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
	{
		faceOffsets += 1 + it->second.facetsNumber();
		for( std::set<UInt>::const_iterator jt = it->second.getFacetsSet().begin(); jt != it->second.getFacetsSet().end(); ++jt )
			faceOffsets += facets.at(*jt).getNumberOfPoints();
		filestr << faceOffsets << std::endl;
	}
	filestr << "\t\t\t\t</DataArray>" << std::endl;
	filestr << "\t\t\t</Cells>" << std::endl;

	filestr << "\t\t</Piece>" << std::endl;
	filestr << "\t</UnstructuredGrid>" << std::endl;
	filestr << "\t</VTKFile>" << std::endl;

	filestr.close();
}

void ExporterVTU::exportFractures(const Mesh3D & mesh, const std::string filename)
{
	std::fstream filestr;

	filestr.open (filename, std::ios_base::out);

	if (filestr.is_open())
	{
		std::cout << std::endl << " File: " << filename << ", successfully opened";
	}
	else
	{
		std::cerr << std::endl << " *** Error: file not opened *** " << std::endl << std::endl;
		return;
	}

	std::cout << std::endl << " Exporting Fractures in Vtu format... " << std::endl;

	const std::vector<Geometry::Point3D> & nodes = mesh.getNodesVector();
	const std::map<UInt,Facet3D> & facets = mesh.getFacetsMap();
	const std::vector<Geometry::Fracture3D> & fractures = mesh.getFn().getNetwork();
	UInt nPoints = nodes.size();
	UInt nCells = 0;
	UInt offsets = 0;

	for(std::vector<Geometry::Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it)
		nCells += it->getNumberOfFractureFacets();

	// Header
	filestr << "<?xml version=\"1.0\"?>" << std::endl;
	filestr << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl;
	filestr << "\t<UnstructuredGrid>" << std::endl;
	filestr << "\t\t<Piece NumberOfPoints=\"" << nPoints << "\" NumberOfCells=\"" << nCells << "\">" << std::endl;

	// CellData
	filestr << "\t\t\t<CellData Scalars=\"cellID\">" << std::endl;
	filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"facetID\" format=\"ascii\">" << std::endl;
	filestr << std::scientific << std::setprecision(0);
	for(std::vector<Geometry::Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
		for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i )
			filestr << it->getFractureFacetsId()[i] << std::endl;
	filestr << "\t\t\t\t</DataArray>" << std::endl;
	filestr << "\t\t\t</CellData>" << std::endl;

	filestr << std::scientific << std::setprecision(10);

	// Points
	filestr << "\t\t\t<Points>" << std::endl;
	filestr << "\t\t\t\t<DataArray type=\"" << "Float32" << "\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
	for( std::vector<Geometry::Point3D>::const_iterator it = nodes.begin(); it != nodes.end(); ++it )
		filestr << it->x() << " " << it->y() << " " << it->z() <<std::endl;
	filestr << "\t\t\t\t</DataArray>" << std::endl;
	filestr << "\t\t\t</Points>" << std::endl;

	// Cells
	filestr << "\t\t\t<Cells>" << std::endl;
	//  Connectivity
	filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
	for(std::vector<Geometry::Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
	{
		for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i)
		{
			for(std::vector<UInt>::const_iterator jt = facets.at(it->getFractureFacetsId()[i]).getVertexesVector().begin(); jt != facets.at(it->getFractureFacetsId()[i]).getVertexesVector().end()-1; ++jt )
				filestr << *jt << " ";
			filestr << *(facets.at(it->getFractureFacetsId()[i]).getVertexesVector().rbegin()) << std::endl;
		}
	}
	filestr << "\t\t\t\t</DataArray>" << std::endl;

	//  Offsets
	filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << std::endl;
	for(std::vector<Geometry::Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
	{
		for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i)
		{
			offsets += facets.at(it->getFractureFacetsId()[i]).getNumberOfPoints();
			filestr << offsets << std::endl;
		}
	}
	filestr << "\t\t\t\t</DataArray>" << std::endl;

	//  Types
	filestr << "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
	for( UInt i=0; i < nCells ; ++i )
		filestr << "7" << std::endl;
	filestr << "\t\t\t\t</DataArray>" << std::endl;

	filestr << "\t\t\t</Cells>" << std::endl;

	filestr << "\t\t</Piece>" << std::endl;
	filestr << "\t</UnstructuredGrid>" << std::endl;
	filestr << "\t</VTKFile>" << std::endl;

	filestr.close();
}

void ExporterVTU::exportMeshWithFractures(const Mesh3D & mesh, const std::string filename)
{
	std::fstream filestr;

	filestr.open (filename.c_str(), std::ios_base::out);

	if (filestr.is_open())
	{
		std::cout << std::endl << " File: " << filename << ", successfully opened";
	}
	else
	{
		std::cerr << std::endl << " *** Error: file not opened *** " << std::endl << std::endl;
		return;
	}

	std::cout << std::endl << " Exporting Mesh3D in Vtu format... " << std::endl;

	const std::vector<Geometry::Point3D> & nodes = mesh.getNodesVector();
	const std::map<UInt, Facet3D> & facets = mesh.getFacetsMap();
	const std::map<UInt, Cell3D> & cells = mesh.getCellsMap();
	const std::vector<Geometry::Fracture3D> & fractures = mesh.getFn().getNetwork();
	UInt nPoints = nodes.size();
	UInt nCells = cells.size();
	UInt nFractures = 0;
	UInt nTotal;
	UInt offsets = 0, faceOffsets = 0;

	for(std::vector<Geometry::Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it)
		nFractures += it->getNumberOfFractureFacets();

	nTotal = nCells + nFractures;

	// Header
	filestr << "<?xml version=\"1.0\"?>" << std::endl;
	filestr << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl;
	filestr << "\t<UnstructuredGrid>" << std::endl;
	filestr << "\t\t<Piece NumberOfPoints=\"" << nPoints << "\" NumberOfCells=\"" << nTotal << "\">" << std::endl;

	// CellData
	filestr << "\t\t\t<CellData Scalars=\"cellID\">" << std::endl;
	filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"cellID\" format=\"ascii\">" << std::endl;
	filestr << std::scientific << std::setprecision(0);
	for( std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
		filestr << it->first << std::endl;
	for(std::vector<Geometry::Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
		for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i )
			filestr << "-1" << std::endl;
	filestr << "\t\t\t\t</DataArray>" << std::endl;

	filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"facetID\" format=\"ascii\">" << std::endl;
	filestr << std::scientific << std::setprecision(0);
	for( std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
		filestr << "-1" << std::endl;
	for(std::vector<Geometry::Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
		for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i )
			filestr << it->getFractureFacetsId()[i] << std::endl;
	filestr << "\t\t\t\t</DataArray>" << std::endl;
	filestr << "\t\t\t</CellData>" << std::endl;

	filestr << std::scientific << std::setprecision(10);

	// Points
	filestr << "\t\t\t<Points>" << std::endl;
	filestr << "\t\t\t\t<DataArray type=\"" << "Float32" << "\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
	for( std::vector<Geometry::Point3D>::const_iterator it = nodes.begin(); it != nodes.end(); ++it )
		filestr << it->x() << " " << it->y() << " " << it->z() <<std::endl;
	filestr << "\t\t\t\t</DataArray>" << std::endl;
	filestr << "\t\t\t</Points>" << std::endl;

	// Cells
	filestr << "\t\t\t<Cells>" << std::endl;
	//  Connectivity
	filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
	for( std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
	{
		for( std::vector<UInt>::const_iterator jt = it->second.getVertexesVector().begin(); jt != it->second.getVertexesVector().end()-1; ++jt )
			filestr << *jt << " ";
		filestr << *(it->second.getVertexesVector().rbegin()) << std::endl;
	}
	for(std::vector<Geometry::Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
	{
		for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i)
		{
			for(std::vector<UInt>::const_iterator jt = facets.at(it->getFractureFacetsId()[i]).getVertexesVector().begin(); jt != facets.at(it->getFractureFacetsId()[i]).getVertexesVector().end()-1; ++jt )
				filestr << *jt << " ";
			filestr << *(facets.at(it->getFractureFacetsId()[i]).getVertexesVector().rbegin()) << std::endl;
		}
	}
	filestr << "\t\t\t\t</DataArray>" << std::endl;

	//  Offsets
	filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << std::endl;
	for( std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
	{
		offsets += it->second.vertexesNumber();
		filestr << offsets << std::endl;
	}
	for(std::vector<Geometry::Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
	{
		for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i)
		{
			offsets += facets.at(it->getFractureFacetsId()[i]).getNumberOfPoints();
			filestr << offsets << std::endl;
		}
	}
	filestr << "\t\t\t\t</DataArray>" << std::endl;

	//  Types
	filestr << "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
	for( UInt i=0; i < nTotal ; ++i )
		filestr << "42" << std::endl;
	filestr << "\t\t\t\t</DataArray>" << std::endl;

	//  Faces
	filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"faces\" format=\"ascii\">" << std::endl;
	for( std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
	{
		filestr << it->second.facetsNumber() << std::endl;
		for( std::set<UInt>::const_iterator jt = it->second.getFacetsSet().begin(); jt != it->second.getFacetsSet().end(); ++jt )
		{
			filestr << facets.at(*jt).getNumberOfPoints() << std::endl;
			for( std::vector<UInt>::const_iterator kt = facets.at(*jt).getVertexesVector().begin(); kt != facets.at(*jt).getVertexesVector().end()-1; ++kt )
				filestr << *kt << " ";
			filestr << *(facets.at(*jt).getVertexesVector().rbegin()) << std::endl;
		}
	}
	for(std::vector<Geometry::Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
	{
		for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i)
		{
			filestr << "1" << std::endl;
			filestr << facets.at(it->getFractureFacetsId()[i]).getNumberOfPoints() << std::endl;
			for(std::vector<UInt>::const_iterator jt = facets.at(it->getFractureFacetsId()[i]).getVertexesVector().begin(); jt != facets.at(it->getFractureFacetsId()[i]).getVertexesVector().end()-1; ++jt )
				filestr << *jt << " ";
			filestr << *(facets.at(it->getFractureFacetsId()[i]).getVertexesVector().rbegin()) << std::endl;
		}
	}
	filestr << "\t\t\t\t</DataArray>" << std::endl;

	//  Faceoffsets
	filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"faceoffsets\" format=\"ascii\">" << std::endl;
	for( std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
	{
		faceOffsets += 1 + it->second.facetsNumber();
		for( std::set<UInt>::const_iterator jt = it->second.getFacetsSet().begin(); jt != it->second.getFacetsSet().end(); ++jt )
			faceOffsets += facets.at(*jt).getNumberOfPoints();
		filestr << faceOffsets << std::endl;
	}
	for(std::vector<Geometry::Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
	{
		for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i)
		{
			faceOffsets += 2 + facets.at(it->getFractureFacetsId()[i]).getNumberOfPoints();
			filestr << faceOffsets << std::endl;
		}
	}
	filestr << "\t\t\t\t</DataArray>" << std::endl;
	filestr << "\t\t\t</Cells>" << std::endl;

	filestr << "\t\t</Piece>" << std::endl;
	filestr << "\t</UnstructuredGrid>" << std::endl;
	filestr << "\t</VTKFile>" << std::endl;

	filestr.close();
}

void ExporterVTU::exportFractureJunctures(const Rigid_Mesh & mesh, const std::string filename)
{
	std::fstream filestr;

	filestr.open (filename.c_str(), std::ios_base::out);

	if (filestr.is_open())
	{
		std::cout << std::endl << " File: " << filename << ", successfully opened";
	}
	else
	{
		std::cerr << std::endl << " *** Error: file not opened *** " << std::endl << std::endl;
		return;
	}

	std::cout << std::endl << " Exporting Fracture Junctures in Vtu format... " << std::endl;

	const std::vector<Geometry::Point3D> & nodes = mesh.getNodesVector();
	const std::vector<Fracture_Facet> & fractures = mesh.getFractureFacetsIdsVector();
	UInt nPoints = 0;
	UInt nCells = 0;
	UInt offsets = 0, faceOffsets = 0;
	UInt count = 0;

	std::set<UInt> localPoints;
	std::set<Fracture_Juncture, Geometry::less< std::pair<UInt,UInt> > > fracturesJunctures;
	std::map<UInt,UInt> GlobalToLocal;

	// Compute # of fracture junctures
	for (auto fractureFacet_it : fractures)
	{
		for(auto juncture_it : fractureFacet_it.getFractureNeighbors())
		{
			localPoints.insert(juncture_it.first.first);
			localPoints.insert(juncture_it.first.second);
			fracturesJunctures.insert( std::make_pair(juncture_it.first.first,juncture_it.first.second) );
		}
	}
	nCells = fracturesJunctures.size();
	nPoints = localPoints.size();

	// Header
	filestr << "<?xml version=\"1.0\"?>" << std::endl;
	filestr << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl;
	filestr << "\t<UnstructuredGrid>" << std::endl;
	filestr << "\t\t<Piece NumberOfPoints=\"" << nPoints << "\" NumberOfCells=\"" << nCells << "\">" << std::endl;

	filestr << std::scientific << std::setprecision(10);

	// Points
	filestr << "\t\t\t<Points>" << std::endl;
	filestr << "\t\t\t\t<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
	for(std::set<UInt>::const_iterator it = localPoints.begin(); it != localPoints.end(); ++it )
	{
		GlobalToLocal.insert( std::make_pair(*it, count++) );
		filestr << nodes[*it].x() << " " << nodes[*it].y() << " " << nodes[*it].z() <<std::endl;
	}
	filestr << "\t\t\t\t</DataArray>" << std::endl;
	filestr << "\t\t\t</Points>" << std::endl;

	// Cells
	filestr << "\t\t\t<Cells>" << std::endl;
	//  Connectivity
	filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
	for(std::set<Fracture_Juncture, Geometry::less< std::pair<UInt,UInt> > >::const_iterator it = fracturesJunctures.begin(); it != fracturesJunctures.end(); ++it )
		filestr << GlobalToLocal[it->first] << " " << GlobalToLocal[it->second] << std::endl;
	filestr << "\t\t\t\t</DataArray>" << std::endl;

	//  Offsets
	filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << std::endl;
	for(std::set<Fracture_Juncture, Geometry::less< std::pair<UInt,UInt> > >::const_iterator it = fracturesJunctures.begin(); it != fracturesJunctures.end(); ++it )
	{
		offsets += 2;
		filestr << offsets << std::endl;
	}
	filestr << "\t\t\t\t</DataArray>" << std::endl;

	//  Types
	filestr << "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
	for( UInt i=0; i < nCells ; ++i )
		filestr << "3" << std::endl;
	filestr << "\t\t\t\t</DataArray>" << std::endl;
	filestr << "\t\t\t</Cells>" << std::endl;

	filestr << "\t\t</Piece>" << std::endl;
	filestr << "\t</UnstructuredGrid>" << std::endl;
	filestr << "\t</VTKFile>" << std::endl;

	filestr.close();
}

void ExporterVTU::exportSolution(const Rigid_Mesh & mesh, const std::string filename, const Eigen::VectorXd & sol)
{
	std::fstream filestr;

	filestr.open (filename.c_str(), std::ios_base::out);

	if (filestr.is_open())
	{
		std::cout << std::endl << " File: " << filename << ", successfully opened";
	}
	else
	{
		std::cerr << std::endl << " *** Error: file not opened *** " << std::endl << std::endl;
		return;
	}

	std::cout << std::endl << " Exporting Solution in Vtu format... " << std::endl;

	const std::vector<Geometry::Point3D> & nodes = mesh.getNodesVector();
	const std::vector<Facet> & facets = mesh.getFacetsVector();
	const std::vector<Cell> & cells = mesh.getCellsVector();
	const std::vector<Fracture_Facet> & fractures = mesh.getFractureFacetsIdsVector();
	UInt nPoints = nodes.size();
	UInt nCells = cells.size();
	UInt nFractures = mesh.getFractureFacetsIdsVector().size();
	UInt nTotal = nCells + nFractures;
	UInt offsets = 0, faceOffsets = 0;

	// Header
	filestr << "<?xml version=\"1.0\"?>" << std::endl;
	filestr << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl;
	filestr << "\t<UnstructuredGrid>" << std::endl;
	filestr << "\t\t<Piece NumberOfPoints=\"" << nPoints << "\" NumberOfCells=\"" << nTotal << "\">" << std::endl;

	// CellData
	filestr << "\t\t\t<CellData Scalars=\"Pressure\">" << std::endl;
	filestr << "\t\t\t\t<DataArray type=\"Float32\" Name=\"Pressure\" format=\"ascii\">" << std::endl;
	filestr << std::scientific << std::setprecision(10);
	for( UInt i = 0; i < sol.size(); ++i )
		filestr << sol[i] << std::endl;
	filestr << "\t\t\t\t</DataArray>" << std::endl;
	filestr << "\t\t\t</CellData>" << std::endl;

	filestr << std::scientific << std::setprecision(10);

	// Points
	filestr << "\t\t\t<Points>" << std::endl;
	filestr << "\t\t\t\t<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
	for( std::vector<Geometry::Point3D>::const_iterator it = nodes.begin(); it != nodes.end(); ++it )
		filestr << it->x() << " " << it->y() << " " << it->z() <<std::endl;
	filestr << "\t\t\t\t</DataArray>" << std::endl;
	filestr << "\t\t\t</Points>" << std::endl;

	// Cells
	filestr << "\t\t\t<Cells>" << std::endl;
	//  Connectivity
	filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
	for( std::vector<Cell>::const_iterator it = cells.begin(); it != cells.end(); ++it )
	{
		for( std::vector<UInt>::const_iterator jt = it->getVertexesIds().begin(); jt != it->getVertexesIds().end()-1; ++jt )
			filestr << *jt << " ";
		filestr << *(it->getVertexesIds().rbegin()) << std::endl;
	}
	for(std::vector<Fracture_Facet>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
	{
		const Facet & f = it->getFacet();
		for(std::vector<UInt>::const_iterator jt = f.getVertexesIds().begin(); jt != f.getVertexesIds().end()-1; ++jt )
			filestr << *jt << " ";
		filestr << *(f.getVertexesIds().rbegin()) << std::endl;

	}
	filestr << "\t\t\t\t</DataArray>" << std::endl;

	//  Offsets
	filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << std::endl;
	for( std::vector<Cell>::const_iterator it = cells.begin(); it != cells.end(); ++it )
	{
		offsets += it->vertexesNumber();
		filestr << offsets << std::endl;
	}
	for(std::vector<Fracture_Facet>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
	{
		offsets += it->getFacet().vertexesNumber();
		filestr << offsets << std::endl;
	}
	filestr << "\t\t\t\t</DataArray>" << std::endl;

	//  Types
	filestr << "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
	for( UInt i=0; i < nTotal ; ++i )
		filestr << "42" << std::endl;
	filestr << "\t\t\t\t</DataArray>" << std::endl;

	//  Faces
	filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"faces\" format=\"ascii\">" << std::endl;
	for( std::vector<Cell>::const_iterator it = cells.begin(); it != cells.end(); ++it )
	{
		filestr << it->facetsNumber() << std::endl;
		for( std::vector<UInt>::const_iterator jt = it->getFacetsIds().begin(); jt != it->getFacetsIds().end(); ++jt )
		{
			filestr << facets[*jt].vertexesNumber() << std::endl;
			for( std::vector<UInt>::const_iterator kt = facets[*jt].getVertexesIds().begin(); kt != facets[*jt].getVertexesIds().end()-1; ++kt )
				filestr << *kt << " ";
			filestr << *(facets[*jt].getVertexesIds().rbegin()) << std::endl;
		}
	}
	for(std::vector<Fracture_Facet>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
	{
		filestr << "1" << std::endl;
		filestr << facets[it->getFacetId()].vertexesNumber() << std::endl;
		for(std::vector<UInt>::const_iterator jt = facets[it->getFacetId()].getVertexesIds().begin(); jt != facets[it->getFacetId()].getVertexesIds().end()-1; ++jt )
			filestr << *jt << " ";
		filestr << *(facets[it->getFacetId()].getVertexesIds().rbegin()) << std::endl;
	}
	filestr << "\t\t\t\t</DataArray>" << std::endl;

	//  Faceoffsets
	filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"faceoffsets\" format=\"ascii\">" << std::endl;
	for( std::vector<Cell>::const_iterator it = cells.begin(); it != cells.end(); ++it )
	{
		faceOffsets += 1 + it->facetsNumber();
		for( std::vector<UInt>::const_iterator jt = it->getFacetsIds().begin(); jt != it->getFacetsIds().end(); ++jt )
			faceOffsets += facets[*jt].vertexesNumber();
		filestr << faceOffsets << std::endl;
	}
	for(std::vector<Fracture_Facet>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
	{
		faceOffsets += 2 + facets[it->getFacetId()].vertexesNumber();
		filestr << faceOffsets << std::endl;
	}
	filestr << "\t\t\t\t</DataArray>" << std::endl;
	filestr << "\t\t\t</Cells>" << std::endl;

	filestr << "\t\t</Piece>" << std::endl;
	filestr << "\t</UnstructuredGrid>" << std::endl;
	filestr << "\t</VTKFile>" << std::endl;

	filestr.close();
}

void ExporterVTU::exportSolutionOnFractures(const Rigid_Mesh & mesh, const std::string filename, const Eigen::VectorXd & sol)
{
	std::fstream filestr;

	filestr.open (filename.c_str(), std::ios_base::out);

	if (filestr.is_open())
	{
		std::cout << std::endl << " File: " << filename << ", successfully opened";
	}
	else
	{
		std::cerr << std::endl << " *** Error: file not opened *** " << std::endl << std::endl;
		return;
	}

	std::cout << std::endl << " Exporting Solution on Fractures in Vtu format... " << std::endl;

	const std::vector<Geometry::Point3D> & nodes = mesh.getNodesVector();
	const std::vector<Facet> & facets = mesh.getFacetsVector();
	const std::vector<Cell> & cells = mesh.getCellsVector();
	const std::vector<Fracture_Facet> & fractures = mesh.getFractureFacetsIdsVector();
	UInt nPoints = nodes.size();
	UInt nCells = cells.size();
	UInt nFractures = mesh.getFractureFacetsIdsVector().size();
	UInt offsets = 0, faceOffsets = 0;

	// Header
	filestr << "<?xml version=\"1.0\"?>" << std::endl;
	filestr << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl;
	filestr << "\t<UnstructuredGrid>" << std::endl;
	filestr << "\t\t<Piece NumberOfPoints=\"" << nPoints << "\" NumberOfCells=\"" << nFractures << "\">" << std::endl;

	// CellData
	filestr << "\t\t\t<CellData Scalars=\"Pressure\">" << std::endl;
	filestr << "\t\t\t\t<DataArray type=\"Float32\" Name=\"Pressure\" format=\"ascii\">" << std::endl;
	filestr << std::scientific << std::setprecision(10);
	for( UInt i = nCells; i < sol.size(); ++i )
		filestr << sol[i] << std::endl;
	filestr << "\t\t\t\t</DataArray>" << std::endl;
	filestr << "\t\t\t</CellData>" << std::endl;

	filestr << std::scientific << std::setprecision(10);

	// Points
	filestr << "\t\t\t<Points>" << std::endl;
	filestr << "\t\t\t\t<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
	for( std::vector<Geometry::Point3D>::const_iterator it = nodes.begin(); it != nodes.end(); ++it )
		filestr << it->x() << " " << it->y() << " " << it->z() <<std::endl;
	filestr << "\t\t\t\t</DataArray>" << std::endl;
	filestr << "\t\t\t</Points>" << std::endl;

	// Cells
	filestr << "\t\t\t<Cells>" << std::endl;
	//  Connectivity
	filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
	for(std::vector<Fracture_Facet>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
	{
		const Facet & f = it->getFacet();
		for(std::vector<UInt>::const_iterator jt = f.getVertexesIds().begin(); jt != f.getVertexesIds().end()-1; ++jt )
			filestr << *jt << " ";
		filestr << *(f.getVertexesIds().rbegin()) << std::endl;

	}
	filestr << "\t\t\t\t</DataArray>" << std::endl;

	//  Offsets
	filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << std::endl;
	for(std::vector<Fracture_Facet>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
	{
		offsets += it->getFacet().vertexesNumber();
		filestr << offsets << std::endl;
	}
	filestr << "\t\t\t\t</DataArray>" << std::endl;

	//  Types
	filestr << "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
	for( UInt i=0; i < nFractures ; ++i )
		filestr << "42" << std::endl;
	filestr << "\t\t\t\t</DataArray>" << std::endl;

	//  Faces
	filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"faces\" format=\"ascii\">" << std::endl;
	for(std::vector<Fracture_Facet>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
	{
		filestr << "1" << std::endl;
		filestr << facets[it->getFacetId()].vertexesNumber() << std::endl;
		for(std::vector<UInt>::const_iterator jt = facets[it->getFacetId()].getVertexesIds().begin(); jt != facets[it->getFacetId()].getVertexesIds().end()-1; ++jt )
			filestr << *jt << " ";
		filestr << *(facets[it->getFacetId()].getVertexesIds().rbegin()) << std::endl;
	}
	filestr << "\t\t\t\t</DataArray>" << std::endl;

	//  Faceoffsets
	filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"faceoffsets\" format=\"ascii\">" << std::endl;
	for(std::vector<Fracture_Facet>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
	{
		faceOffsets += 2 + facets[it->getFacetId()].vertexesNumber();
		filestr << faceOffsets << std::endl;
	}
	filestr << "\t\t\t\t</DataArray>" << std::endl;
	filestr << "\t\t\t</Cells>" << std::endl;

	filestr << "\t\t</Piece>" << std::endl;
	filestr << "\t</UnstructuredGrid>" << std::endl;
	filestr << "\t</VTKFile>" << std::endl;

	filestr.close();
}

void ExporterVTU::exportWithProperties(const Rigid_Mesh & mesh, const std::string filename, const Flag8bit propertiesType, const std::vector<Real> * property)
{
	std::fstream filestr;

	filestr.open (filename.c_str(), std::ios_base::out);

	if (filestr.is_open())
	{
		std::cout << std::endl << " File: " << filename << ", successfully opened";
	}
	else
	{
		std::cerr << std::endl << " *** Error: file not opened *** " << std::endl << std::endl;
		return;
	}

	std::vector<Real> * prop(0);

	std::cout << std::endl << " Exporting Mesh3D in Vtu format... " << std::endl;

	const std::vector<Geometry::Point3D> & nodes = mesh.getNodesVector();
	const std::vector<Facet> & facets = mesh.getFacetsVector();
	const std::vector<Cell> & cells = mesh.getCellsVector();
	const std::vector<Fracture_Facet> & fractures = mesh.getFractureFacetsIdsVector();
	UInt nPoints = nodes.size();
	UInt nCells = cells.size();
	UInt nFractures = mesh.getFractureFacetsIdsVector().size();
	UInt nTotal = nCells + nFractures;
	UInt offsets = 0, faceOffsets = 0;

	if(propertiesType == Other)
		prop = const_cast<std::vector<Real> *>(property);
	if(propertiesType == ZoneCode)
	{
		prop = new std::vector<Real>;
		prop->reserve(nTotal);

		for(std::vector<Cell>::const_iterator it = mesh.getCellsVector().begin() ; it != mesh.getCellsVector().end(); ++it)
			prop->push_back(it->getZoneCode());
		for(std::vector<Fracture_Facet>::const_iterator it = mesh.getFractureFacetsIdsVector().begin() ; it != mesh.getFractureFacetsIdsVector().end(); ++it)
			prop->push_back(it->getZoneCode());
	}
	if(propertiesType == Aperture)
	{
		prop = new std::vector<Real>;
		prop->reserve(nTotal);

		for(std::vector<Cell>::const_iterator it = mesh.getCellsVector().begin() ; it != mesh.getCellsVector().end(); ++it)
			prop->push_back( mesh.getPropertiesMap().getProperties(it->getZoneCode()).M_aperture );
		for(std::vector<Fracture_Facet>::const_iterator it = mesh.getFractureFacetsIdsVector().begin() ; it != mesh.getFractureFacetsIdsVector().end(); ++it)
			prop->push_back( mesh.getPropertiesMap().getProperties(it->getZoneCode()).M_aperture );
	}
	if(propertiesType == Permeability)
	{
		prop = new std::vector<Real>;
		prop->reserve(nTotal);

		for(std::vector<Cell>::const_iterator it = mesh.getCellsVector().begin() ; it != mesh.getCellsVector().end(); ++it)
			prop->push_back( mesh.getPropertiesMap().getProperties(it->getZoneCode()).M_permeability );
		for(std::vector<Fracture_Facet>::const_iterator it = mesh.getFractureFacetsIdsVector().begin() ; it != mesh.getFractureFacetsIdsVector().end(); ++it)
			prop->push_back( mesh.getPropertiesMap().getProperties(it->getZoneCode()).M_permeability );
	}
	if(propertiesType == Porosity)
	{
		prop = new std::vector<Real>;
		prop->reserve(nTotal);

		for(std::vector<Cell>::const_iterator it = mesh.getCellsVector().begin() ; it != mesh.getCellsVector().end(); ++it)
			prop->push_back( mesh.getPropertiesMap().getProperties(it->getZoneCode()).M_porosity );
		for(std::vector<Fracture_Facet>::const_iterator it = mesh.getFractureFacetsIdsVector().begin() ; it != mesh.getFractureFacetsIdsVector().end(); ++it)
			prop->push_back( mesh.getPropertiesMap().getProperties(it->getZoneCode()).M_porosity );
	}
//	if(propertiesType == BorderID)
//	{
//		prop = new std::vector<Real>;
//		prop->reserve(nTotal);
//
//		for(std::vector<Cell>::const_iterator it = mesh.getCellsVector().begin() ; it != mesh.getCellsVector().end(); ++it)
//			prop->push_back( mesh.getPropertiesMap().getProperties(it->getZoneCode()).M_porosity );
//		for(std::vector<Fracture_Facet>::const_iterator it = mesh.getFractureFacetsIdsVector().begin() ; it != mesh.getFractureFacetsIdsVector().end(); ++it)
//			prop->push_back( mesh.getPropertiesMap().getProperties(it->getZoneCode()).M_porosity );
//	}

	// Header
	filestr << "<?xml version=\"1.0\"?>" << std::endl;
	filestr << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl;
	filestr << "\t<UnstructuredGrid>" << std::endl;
	filestr << "\t\t<Piece NumberOfPoints=\"" << nPoints << "\" NumberOfCells=\"" << nTotal << "\">" << std::endl;

	// CellData
	filestr << "\t\t\t<CellData Scalars=\"Prop\">" << std::endl;
	filestr << "\t\t\t\t<DataArray type=\"Float32\" Name=\"Prop\" format=\"ascii\">" << std::endl;
	filestr << std::scientific << std::setprecision(10);
	for( UInt i = 0; i < prop->size(); ++i )
		filestr << (*prop)[i] << std::endl;
	filestr << "\t\t\t\t</DataArray>" << std::endl;
	filestr << "\t\t\t</CellData>" << std::endl;

	filestr << std::scientific << std::setprecision(10);

	// Points
	filestr << "\t\t\t<Points>" << std::endl;
	filestr << "\t\t\t\t<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
	for( std::vector<Geometry::Point3D>::const_iterator it = nodes.begin(); it != nodes.end(); ++it )
		filestr << it->x() << " " << it->y() << " " << it->z() <<std::endl;
	filestr << "\t\t\t\t</DataArray>" << std::endl;
	filestr << "\t\t\t</Points>" << std::endl;

	// Cells
	filestr << "\t\t\t<Cells>" << std::endl;
	//  Connectivity
	filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
	for( std::vector<Cell>::const_iterator it = cells.begin(); it != cells.end(); ++it )
	{
		for( std::vector<UInt>::const_iterator jt = it->getVertexesIds().begin(); jt != it->getVertexesIds().end()-1; ++jt )
			filestr << *jt << " ";
		filestr << *(it->getVertexesIds().rbegin()) << std::endl;
	}
	for(std::vector<Fracture_Facet>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
	{
		const Facet & f = it->getFacet();
		for(std::vector<UInt>::const_iterator jt = f.getVertexesIds().begin(); jt != f.getVertexesIds().end()-1; ++jt )
			filestr << *jt << " ";
		filestr << *(f.getVertexesIds().rbegin()) << std::endl;

	}
	filestr << "\t\t\t\t</DataArray>" << std::endl;

	//  Offsets
	filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << std::endl;
	for( std::vector<Cell>::const_iterator it = cells.begin(); it != cells.end(); ++it )
	{
		offsets += it->vertexesNumber();
		filestr << offsets << std::endl;
	}
	for(std::vector<Fracture_Facet>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
	{
		offsets += it->getFacet().vertexesNumber();
		filestr << offsets << std::endl;
	}
	filestr << "\t\t\t\t</DataArray>" << std::endl;

	//  Types
	filestr << "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
	for( UInt i=0; i < nTotal ; ++i )
		filestr << "42" << std::endl;
	filestr << "\t\t\t\t</DataArray>" << std::endl;

	//  Faces
	filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"faces\" format=\"ascii\">" << std::endl;
	for( std::vector<Cell>::const_iterator it = cells.begin(); it != cells.end(); ++it )
	{
		filestr << it->facetsNumber() << std::endl;
		for( std::vector<UInt>::const_iterator jt = it->getFacetsIds().begin(); jt != it->getFacetsIds().end(); ++jt )
		{
			filestr << facets[*jt].vertexesNumber() << std::endl;
			for( std::vector<UInt>::const_iterator kt = facets[*jt].getVertexesIds().begin(); kt != facets[*jt].getVertexesIds().end()-1; ++kt )
				filestr << *kt << " ";
			filestr << *(facets[*jt].getVertexesIds().rbegin()) << std::endl;
		}
	}
	for(std::vector<Fracture_Facet>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
	{
		filestr << "1" << std::endl;
		filestr << facets[it->getFacetId()].vertexesNumber() << std::endl;
		for(std::vector<UInt>::const_iterator jt = facets[it->getFacetId()].getVertexesIds().begin(); jt != facets[it->getFacetId()].getVertexesIds().end()-1; ++jt )
			filestr << *jt << " ";
		filestr << *(facets[it->getFacetId()].getVertexesIds().rbegin()) << std::endl;
	}
	filestr << "\t\t\t\t</DataArray>" << std::endl;

	//  Faceoffsets
	filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"faceoffsets\" format=\"ascii\">" << std::endl;
	for( std::vector<Cell>::const_iterator it = cells.begin(); it != cells.end(); ++it )
	{
		faceOffsets += 1 + it->facetsNumber();
		for( std::vector<UInt>::const_iterator jt = it->getFacetsIds().begin(); jt != it->getFacetsIds().end(); ++jt )
			faceOffsets += facets[*jt].vertexesNumber();
		filestr << faceOffsets << std::endl;
	}
	for(std::vector<Fracture_Facet>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
	{
		faceOffsets += 2 + facets[it->getFacetId()].vertexesNumber();
		filestr << faceOffsets << std::endl;
	}
	filestr << "\t\t\t\t</DataArray>" << std::endl;
	filestr << "\t\t\t</Cells>" << std::endl;

	filestr << "\t\t</Piece>" << std::endl;
	filestr << "\t</UnstructuredGrid>" << std::endl;
	filestr << "\t</VTKFile>" << std::endl;

	filestr.close();

	if(propertiesType != Other)
		delete prop;

}

void ExporterVTU::exportWithProperties(const Mesh3D & mesh, const std::string filename, const Flag8bit propertiesType, const std::vector<Real> * property)
{
	// TODO
}


void ExporterVTK::exportMesh(const Mesh3D & mesh, const std::string filename)
{
	// TODO
}

void ExporterVTK::exportFractures(const Mesh3D & mesh, const std::string filename)
{
	// TODO
}

void ExporterVTK::exportMeshWithFractures(const Mesh3D & mesh, const std::string filename)
{
	// TODO
}

void ExporterVTK::exportFractureJunctures(const Rigid_Mesh & mesh, const std::string filename)
{
	// TODO
}

void ExporterVTK::exportSolution(const Rigid_Mesh & mesh, const std::string filename, const Eigen::VectorXd & sol)
{
	// TODO
}

void ExporterVTK::exportSolutionOnFractures(const Rigid_Mesh & mesh, const std::string filename, const Eigen::VectorXd & sol)
{
	// TODO
}

void ExporterVTK::exportWithProperties(const Rigid_Mesh & mesh, const std::string filename, const Flag8bit propertiesType, const std::vector<Real> * property)
{
	// TODO
}

void ExporterVTK::exportWithProperties(const Mesh3D & mesh, const std::string filename, const Flag8bit propertiesType, const std::vector<Real> * property)
{
	// TODO
}
