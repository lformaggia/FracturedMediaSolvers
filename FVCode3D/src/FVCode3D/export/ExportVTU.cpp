/*!
 * @file exportVTU.cpp
 * @brief Classes for saving files in the VTU format (definitions).
 */

#include <FVCode3D/export/ExportVTU.hpp>

#include <FVCode3D/mesh/Mesh3D.hpp>
#include <FVCode3D/mesh/RigidMesh.hpp>
#include <FVCode3D/property/Properties.hpp>
#include <FVCode3D/property/Permeability.hpp>
#include <FVCode3D/geometry/Operations.hpp>

namespace FVCode3D
{

void ExporterVTU::exportMesh(const Mesh3D & mesh, const std::string filename) const
{
    std::fstream filestr;

    filestr.open (filename.c_str(), std::ios_base::out);

    if (filestr.is_open())
    {
        std::cout << std::endl << " File: " << filename << ", successfully opened";
    }
    else
    {
        throw std::runtime_error("Error: file " + filename + " not opened.");
    }

    std::cout << std::endl << " Exporting Mesh3D in Vtu format... " << std::endl;

    const std::vector<Point3D> & nodes = mesh.getNodesVector();
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

    filestr << std::scientific << std::setprecision(16);

    // Points
    filestr << "\t\t\t<Points>" << std::endl;
    filestr << "\t\t\t\t<DataArray type=\"" << "Float64" << "\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
    for( std::vector<Point3D>::const_iterator it = nodes.begin(); it != nodes.end(); ++it )
        filestr << it->x() << " " << it->y() << " " << it->z() <<std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;
    filestr << "\t\t\t</Points>" << std::endl;

    // Cells
    filestr << "\t\t\t<Cells>" << std::endl;
    //  Connectivity
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
    for( std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
    {
        for( std::vector<UInt>::const_iterator jt = it->second.getVerticesVector().begin(); jt != it->second.getVerticesVector().end()-1; ++jt )
            filestr << *jt << " ";
        filestr << *(it->second.getVerticesVector().rbegin()) << std::endl;
    }
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    //  Offsets
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << std::endl;
    for( std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
    {
        offsets += it->second.verticesNumber();
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
            filestr << facets.at(*jt).getNumberOfVertices() << std::endl;
            for( std::vector<UInt>::const_iterator kt = facets.at(*jt).getVerticesVector().begin(); kt != facets.at(*jt).getVerticesVector().end()-1; ++kt )
                filestr << *kt << " ";
            filestr << *(facets.at(*jt).getVerticesVector().rbegin()) << std::endl;
        }
    }
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    //  Faceoffsets
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"faceoffsets\" format=\"ascii\">" << std::endl;
    for( std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
    {
        faceOffsets += 1 + it->second.facetsNumber();
        for( std::set<UInt>::const_iterator jt = it->second.getFacetsSet().begin(); jt != it->second.getFacetsSet().end(); ++jt )
            faceOffsets += facets.at(*jt).getNumberOfVertices();
        filestr << faceOffsets << std::endl;
    }
    filestr << "\t\t\t\t</DataArray>" << std::endl;
    filestr << "\t\t\t</Cells>" << std::endl;

    filestr << "\t\t</Piece>" << std::endl;
    filestr << "\t</UnstructuredGrid>" << std::endl;
    filestr << "\t</VTKFile>" << std::endl;

    filestr.close();
}

void ExporterVTU::exportTetrahedralMesh(const Mesh3D & mesh, const std::string filename) const
{
    std::fstream filestr;

    filestr.open (filename.c_str(), std::ios_base::out);

    if (filestr.is_open())
    {
        std::cout << std::endl << " File: " << filename << ", successfully opened";
    }
    else
    {
        throw std::runtime_error("Error: file " + filename + " not opened.");
    }

    std::cout << std::endl << " Exporting Tetrahedral Mesh3D in Vtu format... " << std::endl;

    const std::vector<Point3D> & nodes = mesh.getNodesVector();
    const std::map<UInt, Cell3D> & cells = mesh.getCellsMap();
    UInt nPoints = nodes.size();
    UInt nCells = cells.size();
    UInt offsets = 0;

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

    filestr << std::scientific << std::setprecision(16);

    // Points
    filestr << "\t\t\t<Points>" << std::endl;
    filestr << "\t\t\t\t<DataArray type=\"" << "Float64" << "\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
    for( std::vector<Point3D>::const_iterator it = nodes.begin(); it != nodes.end(); ++it )
        filestr << it->x() << " " << it->y() << " " << it->z() <<std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;
    filestr << "\t\t\t</Points>" << std::endl;

    // Cells
    filestr << "\t\t\t<Cells>" << std::endl;
    //  Connectivity
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
    for( std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
    {
        for( std::vector<UInt>::const_iterator jt = it->second.getVerticesVector().begin(); jt != it->second.getVerticesVector().end()-1; ++jt )
            filestr << *jt << " ";
        filestr << *(it->second.getVerticesVector().rbegin()) << std::endl;
    }
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    //  Offsets
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << std::endl;
    for( std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
    {
        offsets += it->second.verticesNumber();
        filestr << offsets << std::endl;
    }
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    //  Types
    filestr << "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
    for( UInt i=0; i < nCells ; ++i )
        filestr << "10" << std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;
    filestr << "\t\t\t</Cells>" << std::endl;

    filestr << "\t\t</Piece>" << std::endl;
    filestr << "\t</UnstructuredGrid>" << std::endl;
    filestr << "\t</VTKFile>" << std::endl;

    filestr.close();
}

void ExporterVTU::exportFractures(const Mesh3D & mesh, const std::string filename) const
{
    std::fstream filestr;

    filestr.open (filename, std::ios_base::out);

    if (filestr.is_open())
    {
        std::cout << std::endl << " File: " << filename << ", successfully opened";
    }
    else
    {
        throw std::runtime_error("Error: file " + filename + " not opened.");;
    }

    std::cout << std::endl << " Exporting Fractures in Vtu format... " << std::endl;

    const std::vector<Point3D> & nodes = mesh.getNodesVector();
    const std::map<UInt,Facet3D> & facets = mesh.getFacetsMap();
    const std::vector<Fracture3D> & fractures = mesh.getFN().getNetwork();
    UInt nPoints = nodes.size();
    UInt nCells = 0;
    UInt offsets = 0;

    for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it)
        nCells += it->getNumberOfFractureFacets();

    // Header
    filestr << "<?xml version=\"1.0\"?>" << std::endl;
    filestr << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl;
    filestr << "\t<UnstructuredGrid>" << std::endl;
    filestr << "\t\t<Piece NumberOfPoints=\"" << nPoints << "\" NumberOfCells=\"" << nCells << "\">" << std::endl;

    // CellData
    filestr << "\t\t\t<CellData Scalars=\"scalars\">" << std::endl;
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"facetID\" format=\"ascii\">" << std::endl;
    filestr << std::scientific << std::setprecision(0);
    for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
        for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i )
            filestr << it->getFractureFacetsId()[i] << std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"fractureID\" format=\"ascii\">" << std::endl;
    filestr << std::scientific << std::setprecision(0);
    for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
        for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i )
            filestr << it->getId() << std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;
    filestr << "\t\t\t</CellData>" << std::endl;

    filestr << std::scientific << std::setprecision(16);

    // Points
    filestr << "\t\t\t<Points>" << std::endl;
    filestr << "\t\t\t\t<DataArray type=\"" << "Float64" << "\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
    for( std::vector<Point3D>::const_iterator it = nodes.begin(); it != nodes.end(); ++it )
        filestr << it->x() << " " << it->y() << " " << it->z() <<std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;
    filestr << "\t\t\t</Points>" << std::endl;

    // Cells
    filestr << "\t\t\t<Cells>" << std::endl;
    //  Connectivity
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
    for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
    {
        for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i)
        {
            for(std::vector<UInt>::const_iterator jt = facets.at(it->getFractureFacetsId()[i]).getVerticesVector().begin(); jt != facets.at(it->getFractureFacetsId()[i]).getVerticesVector().end()-1; ++jt )
                filestr << *jt << " ";
            filestr << *(facets.at(it->getFractureFacetsId()[i]).getVerticesVector().rbegin()) << std::endl;
        }
    }
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    //  Offsets
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << std::endl;
    for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
    {
        for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i)
        {
            offsets += facets.at(it->getFractureFacetsId()[i]).getNumberOfVertices();
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

void ExporterVTU::exportMeshWithFractures(const Mesh3D & mesh, const std::string filename) const 
{
    std::fstream filestr;

    filestr.open (filename.c_str(), std::ios_base::out);

    if (filestr.is_open())
    {
        std::cout << std::endl << " File: " << filename << ", successfully opened";
    }
    else
    {
        throw std::runtime_error("Error: file " + filename + " not opened.");
    }

    std::cout << std::endl << " Exporting Mesh3D in Vtu format... " << std::endl;

    const std::vector<Point3D> & nodes = mesh.getNodesVector();
    const std::map<UInt, Facet3D> & facets = mesh.getFacetsMap();
    const std::map<UInt, Cell3D> & cells = mesh.getCellsMap();
    const std::vector<Fracture3D> & fractures = mesh.getFN().getNetwork();
    std::set<UInt> ids;
    std::set<UInt>::const_iterator itSet;
    UInt nPoints = nodes.size();
    UInt nCells = cells.size();
    UInt nFractures = 0;
    UInt nTotal;
    UInt offsets = 0, faceOffsets = 0;

    for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it)
        for(UInt i=0; i<it->getNumberOfFractureFacets(); ++i)
          ids.insert(it->getFractureFacetsId()[i]);
    nFractures = ids.size();
    ids.clear();

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
    for(UInt i=0; i<nFractures; ++i )
        filestr << "-1" << std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"facetID\" format=\"ascii\">" << std::endl;
    filestr << std::scientific << std::setprecision(0);
    for( std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
        filestr << "-1" << std::endl;
    for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
    {
        for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i )
        {
            itSet = ids.find(it->getFractureFacetsId()[i]);
            if(itSet == ids.end())
            {
                filestr << it->getFractureFacetsId()[i] << std::endl;
                ids.insert(it->getFractureFacetsId()[i]);
            }
        }
    }
    ids.clear();
    filestr << "\t\t\t\t</DataArray>" << std::endl;
    filestr << "\t\t\t</CellData>" << std::endl;

    filestr << std::scientific << std::setprecision(16);

    // Points
    filestr << "\t\t\t<Points>" << std::endl;
    filestr << "\t\t\t\t<DataArray type=\"" << "Float64" << "\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
    for( std::vector<Point3D>::const_iterator it = nodes.begin(); it != nodes.end(); ++it )
        filestr << it->x() << " " << it->y() << " " << it->z() <<std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;
    filestr << "\t\t\t</Points>" << std::endl;

    // Cells
    filestr << "\t\t\t<Cells>" << std::endl;
    //  Connectivity
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
    for( std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
    {
        for( std::vector<UInt>::const_iterator jt = it->second.getVerticesVector().begin(); jt != it->second.getVerticesVector().end()-1; ++jt )
            filestr << *jt << " ";
        filestr << *(it->second.getVerticesVector().rbegin()) << std::endl;
    }
    for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
    {
        for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i )
        {
            itSet = ids.find(it->getFractureFacetsId()[i]);
            if(itSet == ids.end())
            {
                for(std::vector<UInt>::const_iterator jt = facets.at(it->getFractureFacetsId()[i]).getVerticesVector().begin(); jt != facets.at(it->getFractureFacetsId()[i]).getVerticesVector().end()-1; ++jt )
                    filestr << *jt << " ";
                filestr << *(facets.at(it->getFractureFacetsId()[i]).getVerticesVector().rbegin()) << std::endl;
                ids.insert(it->getFractureFacetsId()[i]);
            }
        }
    }
    ids.clear();
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    //  Offsets
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << std::endl;
    for( std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
    {
        offsets += it->second.verticesNumber();
        filestr << offsets << std::endl;
    }
    for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
    {
        for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i )
        {
            itSet = ids.find(it->getFractureFacetsId()[i]);
            if(itSet == ids.end())
            {
                offsets += facets.at(it->getFractureFacetsId()[i]).getNumberOfVertices();
                filestr << offsets << std::endl;
                ids.insert(it->getFractureFacetsId()[i]);
            }
        }
    }
    ids.clear();
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
            filestr << facets.at(*jt).getNumberOfVertices() << std::endl;
            for( std::vector<UInt>::const_iterator kt = facets.at(*jt).getVerticesVector().begin(); kt != facets.at(*jt).getVerticesVector().end()-1; ++kt )
                filestr << *kt << " ";
            filestr << *(facets.at(*jt).getVerticesVector().rbegin()) << std::endl;
        }
    }
    for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
    {
        for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i )
        {
            itSet = ids.find(it->getFractureFacetsId()[i]);
            if(itSet == ids.end())
            {
                filestr << "1" << std::endl;
                filestr << facets.at(it->getFractureFacetsId()[i]).getNumberOfVertices() << std::endl;
                for(std::vector<UInt>::const_iterator jt = facets.at(it->getFractureFacetsId()[i]).getVerticesVector().begin(); jt != facets.at(it->getFractureFacetsId()[i]).getVerticesVector().end()-1; ++jt )
                    filestr << *jt << " ";
                filestr << *(facets.at(it->getFractureFacetsId()[i]).getVerticesVector().rbegin()) << std::endl;
                ids.insert(it->getFractureFacetsId()[i]);
            }
        }
    }
    ids.clear();
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    //  Faceoffsets
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"faceoffsets\" format=\"ascii\">" << std::endl;
    for( std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
    {
        faceOffsets += 1 + it->second.facetsNumber();
        for( std::set<UInt>::const_iterator jt = it->second.getFacetsSet().begin(); jt != it->second.getFacetsSet().end(); ++jt )
            faceOffsets += facets.at(*jt).getNumberOfVertices();
        filestr << faceOffsets << std::endl;
    }
    for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
    {
        for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i )
        {
            itSet = ids.find(it->getFractureFacetsId()[i]);
            if(itSet == ids.end())
            {
                faceOffsets += 2 + facets.at(it->getFractureFacetsId()[i]).getNumberOfVertices();
                filestr << faceOffsets << std::endl;
                ids.insert(it->getFractureFacetsId()[i]);
            }
        }
    }
    ids.clear();
    filestr << "\t\t\t\t</DataArray>" << std::endl;
    filestr << "\t\t\t</Cells>" << std::endl;

    filestr << "\t\t</Piece>" << std::endl;
    filestr << "\t</UnstructuredGrid>" << std::endl;
    filestr << "\t</VTKFile>" << std::endl;

    filestr.close();
}

void ExporterVTU::exportWireframe(const Mesh3D & mesh, const std::string filename) const
{
    std::fstream filestr;

    filestr.open (filename.c_str(), std::ios_base::out);

    if (filestr.is_open())
    {
        std::cout << std::endl << " File: " << filename << ", successfully opened";
    }
    else
    {
        throw std::runtime_error("Error: file " + filename + " not opened.");
    }

    std::cout << std::endl << " Exporting Wireframe in Vtu format... " << std::endl;

    const std::vector<Point3D> & nodes = mesh.getNodesVector();
    const std::map<UInt, Facet3D> & facets = mesh.getFacetsMap();
    std::map< std::pair<UInt,UInt>, bool > edges;
    std::map< std::pair<UInt,UInt>, bool >::const_iterator itSet;
    UInt edgesFacet = 0;
    UInt nPoints = nodes.size();
    UInt nEdges = 0, offsets = 0;

    for(std::map<UInt, Facet3D>::const_iterator it = facets.begin(); it != facets.end(); ++it)
    {
        edgesFacet = it->second.getNumberOfVertices();
        for(UInt i=0; i<edgesFacet-1; ++i)
        {
            if(it->second.getVertexId(i) < it->second.getVertexId(i+1))
            {
                itSet = edges.find( std::pair<UInt,UInt>(it->second.getVertexId(i), it->second.getVertexId(i+1)) );
                if(itSet == edges.end())
                    edges.insert(std::pair< std::pair<UInt,UInt>, bool>(std::make_pair(it->second.getVertexId(i), it->second.getVertexId(i+1)), it->second.isFracture()));
                else if(it->second.isFracture() == true)
                    edges[std::pair<UInt,UInt>(it->second.getVertexId(i), it->second.getVertexId(i+1))] = true;
            }
            else
            {
                itSet = edges.find( std::pair<UInt,UInt>(it->second.getVertexId(i+1), it->second.getVertexId(i)) );
                if(itSet == edges.end())
                    edges.insert(std::pair< std::pair<UInt,UInt>, bool>(std::make_pair(it->second.getVertexId(i+1), it->second.getVertexId(i)), it->second.isFracture()));
                else if(it->second.isFracture() == true)
                    edges[std::pair<UInt,UInt>(it->second.getVertexId(i+1), it->second.getVertexId(i))] = true;
            }
        }
        if(it->second.getVertexId(edgesFacet-1) < it->second.getVertexId(0))
        {
            itSet = edges.find( std::pair<UInt,UInt>(it->second.getVertexId(edgesFacet-1), it->second.getVertexId(0)) );
            if(itSet == edges.end())
                edges.insert(std::pair< std::pair<UInt,UInt>, bool>(std::make_pair(it->second.getVertexId(edgesFacet-1), it->second.getVertexId(0)), it->second.isFracture()));
            else if(it->second.isFracture() == true)
                edges[std::pair<UInt,UInt>(it->second.getVertexId(edgesFacet-1), it->second.getVertexId(0))] = true;
        }
        else
        {
            itSet = edges.find( std::pair<UInt,UInt>(it->second.getVertexId(0), it->second.getVertexId(edgesFacet-1)) );
            if(itSet == edges.end())
                edges.insert(std::pair< std::pair<UInt,UInt>, bool>(std::make_pair(it->second.getVertexId(0), it->second.getVertexId(edgesFacet-1)), it->second.isFracture()));
            else if(it->second.isFracture() == true)
                edges[std::pair<UInt,UInt>(it->second.getVertexId(0), it->second.getVertexId(edgesFacet-1))] = true;
        }
    }
    nEdges = edges.size();

    // Header
    filestr << "<?xml version=\"1.0\"?>" << std::endl;
    filestr << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl;
    filestr << "\t<UnstructuredGrid>" << std::endl;
    filestr << "\t\t<Piece NumberOfPoints=\"" << nPoints << "\" NumberOfCells=\"" << nEdges << "\">" << std::endl;

    // CellData
    filestr << "\t\t\t<CellData Scalars=\"scalars\">" << std::endl;
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"isFrac\" format=\"ascii\">" << std::endl;
    filestr << std::scientific << std::setprecision(0);
    for(itSet = edges.begin(); itSet != edges.end(); ++itSet )
        filestr << itSet->second << std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;
    filestr << "\t\t\t</CellData>" << std::endl;

    filestr << std::scientific << std::setprecision(16);

    // Points
    filestr << "\t\t\t<Points>" << std::endl;
    filestr << "\t\t\t\t<DataArray type=\"" << "Float64" << "\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
    for( std::vector<Point3D>::const_iterator it = nodes.begin(); it != nodes.end(); ++it )
        filestr << it->x() << " " << it->y() << " " << it->z() <<std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;
    filestr << "\t\t\t</Points>" << std::endl;

    // Cells
    filestr << "\t\t\t<Cells>" << std::endl;
    //  Connectivity
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
    for(itSet = edges.begin(); itSet != edges.end(); ++itSet )
        filestr << itSet->first.first << " " << itSet->first.second << std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    //  Offsets
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << std::endl;
    for(UInt i=0; i < nEdges; ++i )
    {
        offsets += 2;
        filestr << offsets << std::endl;
    }
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    //  Types
    filestr << "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
    for( UInt i=0; i < nEdges ; ++i )
        filestr << "3" << std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    filestr << "\t\t\t</Cells>" << std::endl;

    filestr << "\t\t</Piece>" << std::endl;
    filestr << "\t</UnstructuredGrid>" << std::endl;
    filestr << "\t</VTKFile>" << std::endl;

    filestr.close();
}

void ExporterVTU::exportEdges(const Rigid_Mesh & mesh, const std::string filename) const
{
    std::fstream filestr;

    filestr.open (filename.c_str(), std::ios_base::out);

    if (filestr.is_open())
    {
        std::cout << std::endl << " File: " << filename << ", successfully opened";
    }
    else
    {
        throw std::runtime_error("Error: file " + filename + " not opened.");
    }

    std::cout << std::endl << " Exporting Edges in Vtu format... " << std::endl;

    const std::vector<Point3D> & nodes = mesh.getNodesVector();
    const std::vector<Edge> & edges = mesh.getEdgesVector();
    const std::vector<Regular_Edge> & regularEdges = mesh.getInternalEdgesIdsVector();
    const std::vector<Juncture_Edge> & junctureEdges = mesh.getJunctureEdgesIdsVector();
    const std::vector<Internal_Tip_Edge> & intTipEdges = mesh.getInternalTipEdgesIdsVector();
    const std::vector<Border_Tip_Edge> & borTipEdges = mesh.getBorderTipEdgesIdsVector();
    const std::vector<Pure_Border_Edge> & pureBorEdges = mesh.getPureBorderEdgesIdsVector();

    UInt nPoints = nodes.size();
    UInt nEdges = edges.size();
    UInt offsets = 0;

    // Header
    filestr << "<?xml version=\"1.0\"?>" << std::endl;
    filestr << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl;
    filestr << "\t<UnstructuredGrid>" << std::endl;
    filestr << "\t\t<Piece NumberOfPoints=\"" << nPoints << "\" NumberOfCells=\"" << nEdges << "\">" << std::endl;

    // CellData
    filestr << "\t\t\t<CellData Scalars=\"Prop\">" << std::endl;

    filestr << std::scientific << std::setprecision(0);
    filestr << "\t\t\t\t<DataArray type=\"UInt8\" Name=\"edgeType\" format=\"ascii\">" << std::endl;
    for(UInt i=0; i<regularEdges.size(); ++i)
        filestr << "1" << std::endl;
    for(UInt i=0; i<junctureEdges.size(); ++i)
        filestr << "2" << std::endl;
    for(UInt i=0; i<intTipEdges.size(); ++i)
        filestr << "3" << std::endl;
    for(UInt i=0; i<borTipEdges.size(); ++i)
        filestr << "4" << std::endl;
    for(UInt i=0; i<pureBorEdges.size(); ++i)
        filestr << "5" << std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"edgeId\" format=\"ascii\">" << std::endl;
    for(auto& edge : regularEdges )
        filestr << edge.getId() << std::endl;
    for(auto& edge : junctureEdges )
        filestr << edge.getId() << std::endl;
    for(auto& edge : intTipEdges )
        filestr << edge.getId() << std::endl;
    for(auto& edge : borTipEdges )
        filestr << edge.getId() << std::endl;
    for(auto& edge : pureBorEdges )
        filestr << edge.getId() << std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    filestr << "\t\t\t\t<DataArray type=\"UInt64\" Name=\"edgeBC\" format=\"ascii\">" << std::endl;

    std::for_each( std::begin(regularEdges), std::end(regularEdges),
                   [&filestr] (const Regular_Edge &) { filestr << "0" << std::endl; } );

    std::for_each( std::begin(junctureEdges), std::end(junctureEdges),
                   [&filestr] (const Juncture_Edge &) { filestr << "0" << std::endl; } );

    std::for_each( std::begin(intTipEdges), std::end(intTipEdges),
                   [&filestr] (const Internal_Tip_Edge &) { filestr << "0" << std::endl; } );

    for(auto& edge : borTipEdges )
        filestr << *std::begin(edge.getBorderIds()) << std::endl;
    for(auto& edge : pureBorEdges )
        filestr << *std::begin(edge.getBorderIds()) << std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    filestr << "\t\t\t</CellData>" << std::endl;

    filestr << std::scientific << std::setprecision(16);

    // Points
    filestr << "\t\t\t<Points>" << std::endl;
    filestr << "\t\t\t\t<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
    for( std::vector<Point3D>::const_iterator it = nodes.begin(); it != nodes.end(); ++it )
        filestr << it->x() << " " << it->y() << " " << it->z() <<std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;
    filestr << "\t\t\t</Points>" << std::endl;

    // Cells
    filestr << "\t\t\t<Cells>" << std::endl;
    //  Connectivity
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
    for(auto& edge : regularEdges )
        filestr << edge.getEdge().getEdge().first << " " << edge.getEdge().getEdge().second << std::endl;
    for(auto& edge : junctureEdges )
        filestr << edge.getEdge().getEdge().first << " " << edge.getEdge().getEdge().second << std::endl;
    for(auto& edge : intTipEdges )
        filestr << edge.getEdge().getEdge().first << " " << edge.getEdge().getEdge().second << std::endl;
    for(auto& edge : borTipEdges )
        filestr << edge.getEdge().getEdge().first << " " << edge.getEdge().getEdge().second << std::endl;
    for(auto& edge : pureBorEdges )
        filestr << edge.getEdge().getEdge().first << " " << edge.getEdge().getEdge().second << std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    //  Offsets
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << std::endl;
    for( UInt i=0; i< nEdges; ++i )
    {
        offsets += 2;
        filestr << offsets << std::endl;
    }
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    //  Types
    filestr << "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
    for( UInt i=0; i < nEdges ; ++i )
        filestr << "3" << std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    filestr << "\t\t\t</Cells>" << std::endl;

    filestr << "\t\t</Piece>" << std::endl;
    filestr << "\t</UnstructuredGrid>" << std::endl;
    filestr << "\t</VTKFile>" << std::endl;

    filestr.close();
}

void ExporterVTU::exportFacets(const Rigid_Mesh & mesh, const std::string filename) const
{
    std::fstream filestr;

    filestr.open (filename.c_str(), std::ios_base::out);

    if (filestr.is_open())
    {
        std::cout << std::endl << " File: " << filename << ", successfully opened";
    }
    else
    {
        throw std::runtime_error("Error: file " + filename + " not opened.");
    }

    std::cout << std::endl << " Exporting Facets in Vtu format... " << std::endl;

    const std::vector<Point3D> & nodes = mesh.getNodesVector();
    const std::vector<Facet> & facets = mesh.getFacetsVector();
    UInt nPoints = nodes.size();
    UInt nFacets = facets.size();
    UInt offsets = 0;

    // Header
    filestr << "<?xml version=\"1.0\"?>" << std::endl;
    filestr << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl;
    filestr << "\t<UnstructuredGrid>" << std::endl;
    filestr << "\t\t<Piece NumberOfPoints=\"" << nPoints << "\" NumberOfCells=\"" << nFacets << "\">" << std::endl;

    // CellData
    filestr << "\t\t\t<CellData Scalars=\"facetId\">" << std::endl;

    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"facetId\" format=\"ascii\">" << std::endl;
    filestr << std::scientific << std::setprecision(16);
    for(auto& facet_it : facets)
        filestr << facet_it.getId() << std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"borderId\" format=\"ascii\">" << std::endl;
    filestr << std::scientific << std::setprecision(16);
    for(auto& facet_it : facets)
        filestr << facet_it.getBorderId() << std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    filestr << "\t\t\t</CellData>" << std::endl;

    filestr << std::scientific << std::setprecision(16);

    // Points
    filestr << "\t\t\t<Points>" << std::endl;
    filestr << "\t\t\t\t<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
    for( std::vector<Point3D>::const_iterator it = nodes.begin(); it != nodes.end(); ++it )
        filestr << it->x() << " " << it->y() << " " << it->z() <<std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;
    filestr << "\t\t\t</Points>" << std::endl;

    // Cells
    filestr << "\t\t\t<Cells>" << std::endl;
    //  Connectivity
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
    for(auto& facet_it : facets)
    {
        for(auto vertex_it : facet_it.getVerticesIds())
            filestr << vertex_it << " ";
        filestr << std::endl;
    }
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    //  Offsets
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << std::endl;
    for(auto& facet_it : facets)
    {
        offsets += facet_it.verticesNumber();
        filestr << offsets << std::endl;
    }
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    //  Types
    filestr << "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
    for(UInt i=0; i < nFacets ; ++i)
        filestr << "7" << std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    filestr << "\t\t\t</Cells>" << std::endl;

    filestr << "\t\t</Piece>" << std::endl;
    filestr << "\t</UnstructuredGrid>" << std::endl;
    filestr << "\t</VTKFile>" << std::endl;

    filestr.close();
}

void ExporterVTU::exportFractureJunctures(const Rigid_Mesh & mesh, const std::string filename) const
{
    std::fstream filestr;

    filestr.open (filename.c_str(), std::ios_base::out);

    if (filestr.is_open())
    {
        std::cout << std::endl << " File: " << filename << ", successfully opened";
    }
    else
    {
        throw std::runtime_error("Error: file " + filename + " not opened.");
    }

    std::cout << std::endl << " Exporting Fracture Junctures in Vtu format... " << std::endl;

    const std::vector<Point3D> & nodes = mesh.getNodesVector();
    const std::vector<Fracture_Facet> & fractures = mesh.getFractureFacetsIdsVector();
    UInt nPoints = 0;
    UInt nCells = 0;
    UInt offsets = 0;
    UInt count = 0;

    std::set<UInt> localPoints;
    std::set<Fracture_Juncture, less< std::pair<UInt,UInt> > > fracturesJunctures;
    std::map<UInt,UInt> GlobalToLocal;

    // Compute # of fracture junctures
    for(auto& fractureFacet_it : fractures)
    {
        for(auto& juncture_it : fractureFacet_it.getFractureNeighbors())
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

    filestr << std::scientific << std::setprecision(16);

    // Points
    filestr << "\t\t\t<Points>" << std::endl;
    filestr << "\t\t\t\t<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
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
    for(std::set<Fracture_Juncture, less< std::pair<UInt,UInt> > >::const_iterator it = fracturesJunctures.begin(); it != fracturesJunctures.end(); ++it )
        filestr << GlobalToLocal[it->first] << " " << GlobalToLocal[it->second] << std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    //  Offsets
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << std::endl;
    for(std::set<Fracture_Juncture, less< std::pair<UInt,UInt> > >::const_iterator it = fracturesJunctures.begin(); it != fracturesJunctures.end(); ++it )
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

void ExporterVTU::exportFractureTips(const Rigid_Mesh & mesh, const std::string filename) const
{
    std::fstream filestr;

    filestr.open (filename.c_str(), std::ios_base::out);

    if (filestr.is_open())
    {
        std::cout << std::endl << " File: " << filename << ", successfully opened";
    }
    else
    {
        throw std::runtime_error("Error: file " + filename + " not opened.");
    }

    std::cout << std::endl << " Exporting Fracture Tips in Vtu format... " << std::endl;

    const std::vector<Point3D> & nodes = mesh.getNodesVector();
    const std::vector<Fracture_Facet> & fractures = mesh.getFractureFacetsIdsVector();
    UInt nPoints = 0;
    UInt nCells = 0;
    UInt offsets = 0;
    UInt count = 0;

    std::set<UInt> localPoints;
    std::set<Fracture_Tip, less< std::pair<UInt,UInt> > > fracturesTips;
    std::map<UInt,UInt> GlobalToLocal;

    // Compute # of fracture tips
    for(auto& fractureFacet : fractures)
    {
        for(auto& tips : fractureFacet.getFractureTips())
        {
            localPoints.insert(tips.first);
            localPoints.insert(tips.second);
            fracturesTips.insert( std::make_pair(tips.first, tips.second) );
        }
    }
    nCells = fracturesTips.size();
    nPoints = localPoints.size();

    // Header
    filestr << "<?xml version=\"1.0\"?>" << std::endl;
    filestr << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl;
    filestr << "\t<UnstructuredGrid>" << std::endl;
    filestr << "\t\t<Piece NumberOfPoints=\"" << nPoints << "\" NumberOfCells=\"" << nCells << "\">" << std::endl;

    filestr << std::scientific << std::setprecision(16);

    // Points
    filestr << "\t\t\t<Points>" << std::endl;
    filestr << "\t\t\t\t<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
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
    for(std::set<Fracture_Tip, less< std::pair<UInt,UInt> > >::const_iterator it = fracturesTips.begin(); it != fracturesTips.end(); ++it )
        filestr << GlobalToLocal[it->first] << " " << GlobalToLocal[it->second] << std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    //  Offsets
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << std::endl;
    for(std::set<Fracture_Tip, less< std::pair<UInt,UInt> > >::const_iterator it = fracturesTips.begin(); it != fracturesTips.end(); ++it )
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

void ExporterVTU::exportWithProperties(const Mesh3D & mesh, const PropertiesMap & properties,
    const std::string filename, const Flag16bit propertiesType, const std::vector<Real> * property) const
{
    std::fstream filestr;

    filestr.open (filename.c_str(), std::ios_base::out);

    if (filestr.is_open())
    {
        std::cout << std::endl << " File: " << filename << ", successfully opened";
    }
    else
    {
        throw std::runtime_error("Error: file " + filename + " not opened.");
    }

    std::cout << std::endl << " Exporting Mesh3D in Vtu format... " << std::endl;

    const std::vector<Point3D> & nodes = mesh.getNodesVector();
    const std::map<UInt, Facet3D> & facets = mesh.getFacetsMap();
    const std::map<UInt, Cell3D> & cells = mesh.getCellsMap();
    const std::vector<Fracture3D> & fractures = mesh.getFN().getNetwork();
    std::set<UInt> ids;
    std::set<UInt>::const_iterator itSet;
    UInt nPoints = nodes.size();
    UInt nCells = cells.size();
    UInt nFractures = 0;
    UInt nTotal;
    UInt offsets = 0, faceOffsets = 0;

    for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it)
        for(UInt i=0; i<it->getNumberOfFractureFacets(); ++i)
          ids.insert(it->getFractureFacetsId()[i]);
    nFractures = ids.size();
    ids.clear();

    nTotal = nCells + nFractures;

    // Header
    filestr << "<?xml version=\"1.0\"?>" << std::endl;
    filestr << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl;
    filestr << "\t<UnstructuredGrid>" << std::endl;
    filestr << "\t\t<Piece NumberOfPoints=\"" << nPoints << "\" NumberOfCells=\"" << nTotal << "\">" << std::endl;
    filestr << std::scientific << std::setprecision(16);

    // CellData
    filestr << "\t\t\t<CellData Scalars=\"Prop\">" << std::endl;
    if(propertiesType & ZoneCode)
    {
        filestr << std::scientific << std::setprecision(0);
        filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"zoneCode\" format=\"ascii\">" << std::endl;
        for(std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
            filestr << it->second.getZoneCode() << std::endl;
        for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
        {
            for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i )
            {
                itSet = ids.find(it->getFractureFacetsId()[i]);
                if(itSet == ids.end())
                {
                    filestr << facets.at(it->getFractureFacetsId()[i]).getZoneCode() << std::endl;
                    ids.insert(it->getFractureFacetsId()[i]);
                }
            }
        }
        ids.clear();
        filestr << "\t\t\t\t</DataArray>" << std::endl;
    }
    if(propertiesType & IsFrac)
    {
        filestr << std::scientific << std::setprecision(0);
        filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"isFrac\" format=\"ascii\">" << std::endl;
        for(UInt i=0; i < nCells; ++i )
            filestr << "0" << std::endl;
        for(UInt i=0; i < nFractures; ++i )
            filestr << "1" << std::endl;
        filestr << "\t\t\t\t</DataArray>" << std::endl;
    }
    if(propertiesType & Aperture)
    {
        filestr << std::scientific << std::setprecision(16);
        filestr << "\t\t\t\t<DataArray type=\"Float64\" Name=\"aperture\" format=\"ascii\">" << std::endl;
        for(UInt i=0; i < nCells; ++i )
            filestr << "-1" << std::endl;
        for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
        {
            for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i )
            {
                itSet = ids.find(it->getFractureFacetsId()[i]);
                if(itSet == ids.end())
                {
                    filestr << properties.getProperties(facets.at(it->getFractureFacetsId()[i]).getZoneCode()).M_aperture << std::endl;
                    ids.insert(it->getFractureFacetsId()[i]);
                }
            }
        }
        ids.clear();
        filestr << "\t\t\t\t</DataArray>" << std::endl;
    }
    if(propertiesType & Porosity)
    {
        filestr << std::scientific << std::setprecision(16);
        filestr << "\t\t\t\t<DataArray type=\"Float64\" Name=\"porosity\" format=\"ascii\">" << std::endl;
        for(std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
            filestr << properties.getProperties(it->second.getZoneCode()).M_porosity << std::endl;
        for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
        {
            for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i )
            {
                itSet = ids.find(it->getFractureFacetsId()[i]);
                if(itSet == ids.end())
                {
                    filestr << properties.getProperties(facets.at(it->getFractureFacetsId()[i]).getZoneCode()).M_porosity << std::endl;
                    ids.insert(it->getFractureFacetsId()[i]);
                }
            }
        }
        ids.clear();
        filestr << "\t\t\t\t</DataArray>" << std::endl;
    }
    if(propertiesType & Permeability)
    {
        filestr << std::scientific << std::setprecision(16);
        filestr << "\t\t\t\t<DataArray type=\"Float64\" Name=\"permeability\" format=\"ascii\">" << std::endl;
        for(std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
            filestr << properties.getProperties(it->second.getZoneCode()).M_permeability->norm() << std::endl;
        for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
        {
            for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i )
            {
                itSet = ids.find(it->getFractureFacetsId()[i]);
                if(itSet == ids.end())
                {
                      filestr <<
                      properties.getProperties(facets.at(it->getFractureFacetsId()[i]).getZoneCode()).M_permeability->norm() << std::endl;
                      ids.insert(it->getFractureFacetsId()[i]);
                }
            }
        }
        ids.clear();
        filestr << "\t\t\t\t</DataArray>" << std::endl;
    }
    if(propertiesType & BorderID)
    {
        filestr << std::scientific << std::setprecision(0);
        filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"borderID\" format=\"ascii\">" << std::endl;
        for(UInt i=0; i < nCells; ++i )
            filestr << "0" << std::endl;
        for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
        {
            for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i )
            {
                itSet = ids.find(it->getFractureFacetsId()[i]);
                if(itSet == ids.end())
                {
                    filestr << facets.at(it->getFractureFacetsId()[i]).getBorderId() << std::endl;
                    ids.insert(it->getFractureFacetsId()[i]);
                }
            }
        }
        ids.clear();
        filestr << "\t\t\t\t</DataArray>" << std::endl;
    }
    if(propertiesType & ElementID)
    {
        filestr << std::scientific << std::setprecision(0);
        filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"cellID\" format=\"ascii\">" << std::endl;
        for(std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
            filestr << it->first << std::endl;
        for(UInt i=0; i < nFractures; ++i )
            filestr << "-1" << std::endl;
        filestr << "\t\t\t\t</DataArray>" << std::endl;
        filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"facetID\" format=\"ascii\">" << std::endl;
        for(UInt i=0; i < nCells; ++i )
            filestr << "-1" << std::endl;
        for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
        {
            for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i )
            {
                itSet = ids.find(it->getFractureFacetsId()[i]);
                if(itSet == ids.end())
                {
                    filestr << it->getFractureFacetsId()[i] << std::endl;
                    ids.insert(it->getFractureFacetsId()[i]);
                }
            }
        }
        ids.clear();
        filestr << "\t\t\t\t</DataArray>" << std::endl;
    }
    if(propertiesType & Other)
    {
        filestr << std::scientific << std::setprecision(16);
        filestr << "\t\t\t\t<DataArray type=\"Float64\" Name=\"prop\" format=\"ascii\">" << std::endl;
        for( UInt i=0; i < nTotal ; ++i )
            filestr << (*property)[i] << std::endl;
        filestr << "\t\t\t\t</DataArray>" << std::endl;
    }
    filestr << "\t\t\t</CellData>" << std::endl;

    filestr << std::scientific << std::setprecision(16);

    // Points
    filestr << "\t\t\t<Points>" << std::endl;
    filestr << "\t\t\t\t<DataArray type=\"" << "Float64" << "\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
    for( std::vector<Point3D>::const_iterator it = nodes.begin(); it != nodes.end(); ++it )
        filestr << it->x() << " " << it->y() << " " << it->z() <<std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;
    filestr << "\t\t\t</Points>" << std::endl;

    // Cells
    filestr << "\t\t\t<Cells>" << std::endl;
    //  Connectivity
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
    for( std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
    {
        for( std::vector<UInt>::const_iterator jt = it->second.getVerticesVector().begin(); jt != it->second.getVerticesVector().end()-1; ++jt )
            filestr << *jt << " ";
        filestr << *(it->second.getVerticesVector().rbegin()) << std::endl;
    }
    for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
    {
        for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i)
        {
            itSet = ids.find(it->getFractureFacetsId()[i]);
            if(itSet == ids.end())
            {
                for(std::vector<UInt>::const_iterator jt = facets.at(it->getFractureFacetsId()[i]).getVerticesVector().begin(); jt != facets.at(it->getFractureFacetsId()[i]).getVerticesVector().end()-1; ++jt )
                    filestr << *jt << " ";
                filestr << *(facets.at(it->getFractureFacetsId()[i]).getVerticesVector().rbegin()) << std::endl;
                ids.insert(it->getFractureFacetsId()[i]);
            }
        }
    }
    ids.clear();
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    //  Offsets
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << std::endl;
    for( std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
    {
        offsets += it->second.verticesNumber();
        filestr << offsets << std::endl;
    }
    for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
    {
        for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i)
        {
            itSet = ids.find(it->getFractureFacetsId()[i]);
            if(itSet == ids.end())
            {
                offsets += facets.at(it->getFractureFacetsId()[i]).getNumberOfVertices();
                filestr << offsets << std::endl;
                ids.insert(it->getFractureFacetsId()[i]);
            }
        }
    }
    ids.clear();
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
            filestr << facets.at(*jt).getNumberOfVertices() << std::endl;
            for( std::vector<UInt>::const_iterator kt = facets.at(*jt).getVerticesVector().begin(); kt != facets.at(*jt).getVerticesVector().end()-1; ++kt )
                filestr << *kt << " ";
            filestr << *(facets.at(*jt).getVerticesVector().rbegin()) << std::endl;
        }
    }
    for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
    {
        for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i)
        {
            itSet = ids.find(it->getFractureFacetsId()[i]);
            if(itSet == ids.end())
            {
                filestr << "1" << std::endl;
                filestr << facets.at(it->getFractureFacetsId()[i]).getNumberOfVertices() << std::endl;
                for(std::vector<UInt>::const_iterator jt = facets.at(it->getFractureFacetsId()[i]).getVerticesVector().begin(); jt != facets.at(it->getFractureFacetsId()[i]).getVerticesVector().end()-1; ++jt )
                    filestr << *jt << " ";
                filestr << *(facets.at(it->getFractureFacetsId()[i]).getVerticesVector().rbegin()) << std::endl;
                ids.insert(it->getFractureFacetsId()[i]);
            }
        }
    }
    ids.clear();
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    //  Faceoffsets
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"faceoffsets\" format=\"ascii\">" << std::endl;
    for( std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
    {
        faceOffsets += 1 + it->second.facetsNumber();
        for( std::set<UInt>::const_iterator jt = it->second.getFacetsSet().begin(); jt != it->second.getFacetsSet().end(); ++jt )
            faceOffsets += facets.at(*jt).getNumberOfVertices();
        filestr << faceOffsets << std::endl;
    }
    for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
    {
        for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i)
        {
            itSet = ids.find(it->getFractureFacetsId()[i]);
            if(itSet == ids.end())
            {
                faceOffsets += 2 + facets.at(it->getFractureFacetsId()[i]).getNumberOfVertices();
                filestr << faceOffsets << std::endl;
                ids.insert(it->getFractureFacetsId()[i]);
            }
        }
    }
    ids.clear();
    filestr << "\t\t\t\t</DataArray>" << std::endl;
    filestr << "\t\t\t</Cells>" << std::endl;

    filestr << "\t\t</Piece>" << std::endl;
    filestr << "\t</UnstructuredGrid>" << std::endl;
    filestr << "\t</VTKFile>" << std::endl;

    filestr.close();
}

void ExporterVTU::exportWithProperties(const Mesh3D & mesh, const PropertiesMap & properties, const std::string filename) const
{
    std::fstream filestr;

    filestr.open (filename.c_str(), std::ios_base::out);

    if (filestr.is_open())
    {
        std::cout << std::endl << " File: " << filename << ", successfully opened";
    }
    else
    {
        throw std::runtime_error("Error: file " + filename + " not opened.");
    }

    std::cout << std::endl << " Exporting Mesh3D in Vtu format... " << std::endl;

    const std::vector<Point3D> & nodes = mesh.getNodesVector();
    const std::map<UInt, Facet3D> & facets = mesh.getFacetsMap();
    const std::map<UInt, Cell3D> & cells = mesh.getCellsMap();
    const std::vector<Fracture3D> & fractures = mesh.getFN().getNetwork();
    std::set<UInt> ids;
    std::set<UInt>::const_iterator itSet;
    UInt nPoints = nodes.size();
    UInt nCells = cells.size();
    UInt nFractures = 0;
    UInt nTotal;
    UInt offsets = 0, faceOffsets = 0;

    for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it)
        for(UInt i=0; i<it->getNumberOfFractureFacets(); ++i)
          ids.insert(it->getFractureFacetsId()[i]);
    nFractures = ids.size();
    ids.clear();

    nTotal = nCells + nFractures;

    // Header
    filestr << "<?xml version=\"1.0\"?>" << std::endl;
    filestr << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl;
    filestr << "\t<UnstructuredGrid>" << std::endl;
    filestr << "\t\t<Piece NumberOfPoints=\"" << nPoints << "\" NumberOfCells=\"" << nTotal << "\">" << std::endl;
    filestr << std::scientific << std::setprecision(16);

    // CellData
    filestr << "\t\t\t<CellData Scalars=\"PropScal\">" << std::endl;
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"isFrac\" format=\"ascii\">" << std::endl;
    for(UInt i=0; i < nCells; ++i )
        filestr << "0" << std::endl;
    for(UInt i=0; i < nFractures; ++i )
        filestr << "1" << std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    filestr << "\t\t\t\t<DataArray type=\"Float64\" Name=\"aperture\" format=\"ascii\">" << std::endl;
    for(UInt i=0; i < nCells; ++i )
        filestr << "-1" << std::endl;
    for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
    {
        for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i )
        {
            itSet = ids.find(it->getFractureFacetsId()[i]);
            if(itSet == ids.end())
            {
                filestr << properties.getProperties(facets.at(it->getFractureFacetsId()[i]).getZoneCode()).M_aperture << std::endl;
                ids.insert(it->getFractureFacetsId()[i]);
            }
        }
    }
    ids.clear();
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    filestr << "\t\t\t\t<DataArray type=\"Float64\" Name=\"porosity\" format=\"ascii\">" << std::endl;
    for(std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
        filestr << properties.getProperties(it->second.getZoneCode()).M_porosity << std::endl;
    for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
    {
        for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i )
        {
            itSet = ids.find(it->getFractureFacetsId()[i]);
            if(itSet == ids.end())
            {
                  filestr << properties.getProperties(facets.at(it->getFractureFacetsId()[i]).getZoneCode()).M_porosity << std::endl;
                  ids.insert(it->getFractureFacetsId()[i]);
            }
        }
    }
    ids.clear();
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    filestr << "\t\t\t\t<DataArray type=\"Float64\" Name=\"permeability\" format=\"ascii\">" << std::endl;
    for(std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
        filestr << properties.getProperties(it->second.getZoneCode()).M_permeability->norm() << std::endl;
    for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
    {
        for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i )
        {
            itSet = ids.find(it->getFractureFacetsId()[i]);
            if(itSet == ids.end())
            {
                  filestr <<
                  properties.getProperties(facets.at(it->getFractureFacetsId()[i]).getZoneCode()).M_permeability->norm() << std::endl;
                  ids.insert(it->getFractureFacetsId()[i]);
            }
        }
    }
    ids.clear();
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    filestr << std::scientific << std::setprecision(0);

    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"cellID\" format=\"ascii\">" << std::endl;
    for(std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
        filestr << it->first << std::endl;
    for(UInt i=0; i < nFractures; ++i )
        filestr << "-1" << std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"facetID\" format=\"ascii\">" << std::endl;
    for(UInt i=0; i < nCells; ++i )
        filestr << "-1" << std::endl;
    for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
    {
        for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i )
        {
            itSet = ids.find(it->getFractureFacetsId()[i]);
            if(itSet == ids.end())
            {
                filestr << it->getFractureFacetsId()[i] << std::endl;
                ids.insert(it->getFractureFacetsId()[i]);
            }
        }
    }
    ids.clear();
    filestr << "\t\t\t\t</DataArray>" << std::endl;
    filestr << "\t\t\t</CellData>" << std::endl;

    filestr << std::scientific << std::setprecision(16);

    // Points
    filestr << "\t\t\t<Points>" << std::endl;
    filestr << "\t\t\t\t<DataArray type=\"" << "Float64" << "\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
    for( std::vector<Point3D>::const_iterator it = nodes.begin(); it != nodes.end(); ++it )
        filestr << it->x() << " " << it->y() << " " << it->z() <<std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;
    filestr << "\t\t\t</Points>" << std::endl;

    // Cells
    filestr << "\t\t\t<Cells>" << std::endl;
    //  Connectivity
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
    for( std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
    {
        for( std::vector<UInt>::const_iterator jt = it->second.getVerticesVector().begin(); jt != it->second.getVerticesVector().end()-1; ++jt )
            filestr << *jt << " ";
        filestr << *(it->second.getVerticesVector().rbegin()) << std::endl;
    }
    for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
    {
        for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i)
        {
            itSet = ids.find(it->getFractureFacetsId()[i]);
            if(itSet == ids.end())
            {
                for(std::vector<UInt>::const_iterator jt = facets.at(it->getFractureFacetsId()[i]).getVerticesVector().begin(); jt != facets.at(it->getFractureFacetsId()[i]).getVerticesVector().end()-1; ++jt )
                    filestr << *jt << " ";
                filestr << *(facets.at(it->getFractureFacetsId()[i]).getVerticesVector().rbegin()) << std::endl;
                ids.insert(it->getFractureFacetsId()[i]);
            }
        }
    }
    ids.clear();
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    //  Offsets
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << std::endl;
    for( std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
    {
        offsets += it->second.verticesNumber();
        filestr << offsets << std::endl;
    }
    for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
    {
        for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i)
        {
            itSet = ids.find(it->getFractureFacetsId()[i]);
            if(itSet == ids.end())
            {
                offsets += facets.at(it->getFractureFacetsId()[i]).getNumberOfVertices();
                filestr << offsets << std::endl;
                ids.insert(it->getFractureFacetsId()[i]);
            }
        }
    }
    ids.clear();
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
            filestr << facets.at(*jt).getNumberOfVertices() << std::endl;
            for( std::vector<UInt>::const_iterator kt = facets.at(*jt).getVerticesVector().begin(); kt != facets.at(*jt).getVerticesVector().end()-1; ++kt )
                filestr << *kt << " ";
            filestr << *(facets.at(*jt).getVerticesVector().rbegin()) << std::endl;
        }
    }
    for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
    {
        for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i)
        {
            itSet = ids.find(it->getFractureFacetsId()[i]);
            if(itSet == ids.end())
            {
                filestr << "1" << std::endl;
                filestr << facets.at(it->getFractureFacetsId()[i]).getNumberOfVertices() << std::endl;
                for(std::vector<UInt>::const_iterator jt = facets.at(it->getFractureFacetsId()[i]).getVerticesVector().begin(); jt != facets.at(it->getFractureFacetsId()[i]).getVerticesVector().end()-1; ++jt )
                    filestr << *jt << " ";
                filestr << *(facets.at(it->getFractureFacetsId()[i]).getVerticesVector().rbegin()) << std::endl;
                ids.insert(it->getFractureFacetsId()[i]);
            }
        }
    }
    ids.clear();
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    //  Faceoffsets
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"faceoffsets\" format=\"ascii\">" << std::endl;
    for( std::map<UInt,Cell3D>::const_iterator it = cells.begin(); it != cells.end(); ++it )
    {
        faceOffsets += 1 + it->second.facetsNumber();
        for( std::set<UInt>::const_iterator jt = it->second.getFacetsSet().begin(); jt != it->second.getFacetsSet().end(); ++jt )
            faceOffsets += facets.at(*jt).getNumberOfVertices();
        filestr << faceOffsets << std::endl;
    }
    for(std::vector<Fracture3D>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
    {
        for(UInt i=0; i < it->getNumberOfFractureFacets(); ++i)
        {
            itSet = ids.find(it->getFractureFacetsId()[i]);
            if(itSet == ids.end())
            {
                faceOffsets += 2 + facets.at(it->getFractureFacetsId()[i]).getNumberOfVertices();
                filestr << faceOffsets << std::endl;
                ids.insert(it->getFractureFacetsId()[i]);
            }
        }
    }
    ids.clear();
    filestr << "\t\t\t\t</DataArray>" << std::endl;
    filestr << "\t\t\t</Cells>" << std::endl;

    filestr << "\t\t</Piece>" << std::endl;
    filestr << "\t</UnstructuredGrid>" << std::endl;
    filestr << "\t</VTKFile>" << std::endl;

    filestr.close();
}

void ExporterVTU::exportWithProperties(const Rigid_Mesh & mesh, const std::string filename, const Flag16bit propertiesType, const std::vector<Real> * property) const
{
    std::fstream filestr;

    filestr.open (filename.c_str(), std::ios_base::out);

    if (filestr.is_open())
    {
        std::cout << std::endl << " File: " << filename << ", successfully opened";
    }
    else
    {
        throw std::runtime_error("Error: file " + filename + " not opened.");
    }

    std::cout << std::endl << " Exporting Mesh3D in Vtu format... " << std::endl;

    const std::vector<Point3D> & nodes = mesh.getNodesVector();
    const std::vector<Facet> & facets = mesh.getFacetsVector();
    const std::vector<Cell> & cells = mesh.getCellsVector();
    const std::vector<Fracture_Facet> & fractures = mesh.getFractureFacetsIdsVector();
    const PropertiesMap & properties = mesh.getPropertiesMap();
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
    filestr << "\t\t\t<CellData Scalars=\"Prop\">" << std::endl;
    if(propertiesType & ZoneCode)
    {
        filestr << std::scientific << std::setprecision(0);
        filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"zoneCode\" format=\"ascii\">" << std::endl;
        for(std::vector<Cell>::const_iterator it = cells.begin(); it != cells.end(); ++it )
            filestr << it->getZoneCode() << std::endl;
        for(std::vector<Fracture_Facet>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
            filestr << it->getZoneCode() << std::endl;
        filestr << "\t\t\t\t</DataArray>" << std::endl;
    }
    if(propertiesType & IsFrac)
    {
        filestr << std::scientific << std::setprecision(0);
        filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"isFrac\" format=\"ascii\">" << std::endl;
        for(UInt i=0; i < nCells; ++i )
            filestr << "0" << std::endl;
        for(UInt i=0; i < nFractures; ++i )
            filestr << "1" << std::endl;
        filestr << "\t\t\t\t</DataArray>" << std::endl;
    }
    if(propertiesType & Aperture)
    {
        filestr << std::scientific << std::setprecision(16);
        filestr << "\t\t\t\t<DataArray type=\"Float64\" Name=\"aperture\" format=\"ascii\">" << std::endl;
        for(UInt i=0; i < nCells; ++i )
            filestr << "-1" << std::endl;
        for(std::vector<Fracture_Facet>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
            filestr << properties.getProperties(it->getZoneCode()).M_aperture << std::endl;
        filestr << "\t\t\t\t</DataArray>" << std::endl;
    }
    if(propertiesType & Porosity)
    {
        filestr << std::scientific << std::setprecision(16);
        filestr << "\t\t\t\t<DataArray type=\"Float64\" Name=\"porosity\" format=\"ascii\">" << std::endl;
        for(std::vector<Cell>::const_iterator it = cells.begin(); it != cells.end(); ++it )
            filestr << properties.getProperties(it->getZoneCode()).M_porosity << std::endl;
        for(std::vector<Fracture_Facet>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
            filestr <<  properties.getProperties(it->getZoneCode()).M_porosity << std::endl;
        filestr << "\t\t\t\t</DataArray>" << std::endl;
    }
    if(propertiesType & Permeability)
    {
        filestr << std::scientific << std::setprecision(16);
        filestr << "\t\t\t\t<DataArray type=\"Float64\" Name=\"permeability\" format=\"ascii\">" << std::endl;
        for(std::vector<Cell>::const_iterator it = cells.begin(); it != cells.end(); ++it )
            filestr << properties.getProperties(it->getZoneCode()).M_permeability->norm() << std::endl;
        for(std::vector<Fracture_Facet>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
            filestr <<  properties.getProperties(it->getZoneCode()).M_permeability->norm() << std::endl;
        filestr << "\t\t\t\t</DataArray>" << std::endl;
    }
    if(propertiesType & BorderID)
    {
        filestr << std::scientific << std::setprecision(0);
        filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"borderID\" format=\"ascii\">" << std::endl;
        for(UInt i=0; i < nCells; ++i )
            filestr << "0" << std::endl;
        for(UInt i=0; i < nFractures; ++i )
            filestr <<  "0" << std::endl;
        filestr << "\t\t\t\t</DataArray>" << std::endl;
    }
    if(propertiesType & ElementID)
    {
        filestr << std::scientific << std::setprecision(0);
        filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"cellID\" format=\"ascii\">" << std::endl;
        for(std::vector<Cell>::const_iterator it = cells.begin(); it != cells.end(); ++it )
            filestr << it->getId() << std::endl;
        for(UInt i=0; i < nFractures; ++i )
            filestr << "-1" << std::endl;
        filestr << "\t\t\t\t</DataArray>" << std::endl;
        filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"facetID\" format=\"ascii\">" << std::endl;
        for(UInt i=0; i < nCells; ++i )
            filestr << "-1" << std::endl;
        for(std::vector<Fracture_Facet>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
        {
            filestr << it->getFractureId() << std::endl;
        }
        filestr << "\t\t\t\t</DataArray>" << std::endl;
    }
    if(propertiesType & Other)
    {
        filestr << std::scientific << std::setprecision(16);
        filestr << "\t\t\t\t<DataArray type=\"Float64\" Name=\"prop\" format=\"ascii\">" << std::endl;
        for( UInt i=0; i < nTotal ; ++i )
            filestr << (*property)[i] << std::endl;
        filestr << "\t\t\t\t</DataArray>" << std::endl;
    }
    filestr << "\t\t\t</CellData>" << std::endl;

    filestr << std::scientific << std::setprecision(16);

    // Points
    filestr << "\t\t\t<Points>" << std::endl;
    filestr << "\t\t\t\t<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
    for( std::vector<Point3D>::const_iterator it = nodes.begin(); it != nodes.end(); ++it )
        filestr << it->x() << " " << it->y() << " " << it->z() <<std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;
    filestr << "\t\t\t</Points>" << std::endl;

    // Cells
    filestr << "\t\t\t<Cells>" << std::endl;
    //  Connectivity
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
    for( std::vector<Cell>::const_iterator it = cells.begin(); it != cells.end(); ++it )
    {
        for( std::vector<UInt>::const_iterator jt = it->getVerticesIds().begin(); jt != it->getVerticesIds().end()-1; ++jt )
            filestr << *jt << " ";
        filestr << *(it->getVerticesIds().rbegin()) << std::endl;
    }
    for(std::vector<Fracture_Facet>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
    {
        const Facet & f = it->getFacet();
        for(std::vector<UInt>::const_iterator jt = f.getVerticesIds().begin(); jt != f.getVerticesIds().end()-1; ++jt )
            filestr << *jt << " ";
        filestr << *(f.getVerticesIds().rbegin()) << std::endl;

    }
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    //  Offsets
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << std::endl;
    for( std::vector<Cell>::const_iterator it = cells.begin(); it != cells.end(); ++it )
    {
        offsets += it->verticesNumber();
        filestr << offsets << std::endl;
    }
    for(std::vector<Fracture_Facet>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
    {
        offsets += it->getFacet().verticesNumber();
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
            filestr << facets[*jt].verticesNumber() << std::endl;
            for( std::vector<UInt>::const_iterator kt = facets[*jt].getVerticesIds().begin(); kt != facets[*jt].getVerticesIds().end()-1; ++kt )
                filestr << *kt << " ";
            filestr << *(facets[*jt].getVerticesIds().rbegin()) << std::endl;
        }
    }
    for(std::vector<Fracture_Facet>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
    {
        filestr << "1" << std::endl;
        filestr << facets[it->getId()].verticesNumber() << std::endl;
        for(std::vector<UInt>::const_iterator jt = facets[it->getId()].getVerticesIds().begin(); jt != facets[it->getId()].getVerticesIds().end()-1; ++jt )
            filestr << *jt << " ";
        filestr << *(facets[it->getId()].getVerticesIds().rbegin()) << std::endl;
    }
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    //  Faceoffsets
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"faceoffsets\" format=\"ascii\">" << std::endl;
    for( std::vector<Cell>::const_iterator it = cells.begin(); it != cells.end(); ++it )
    {
        faceOffsets += 1 + it->facetsNumber();
        for( std::vector<UInt>::const_iterator jt = it->getFacetsIds().begin(); jt != it->getFacetsIds().end(); ++jt )
            faceOffsets += facets[*jt].verticesNumber();
        filestr << faceOffsets << std::endl;
    }
    for(std::vector<Fracture_Facet>::const_iterator it = fractures.begin(); it != fractures.end(); ++it )
    {
        faceOffsets += 2 + facets[it->getId()].verticesNumber();
        filestr << faceOffsets << std::endl;
    }
    filestr << "\t\t\t\t</DataArray>" << std::endl;
    filestr << "\t\t\t</Cells>" << std::endl;

    filestr << "\t\t</Piece>" << std::endl;
    filestr << "\t</UnstructuredGrid>" << std::endl;
    filestr << "\t</VTKFile>" << std::endl;

    filestr.close();
}

} // namespace FVCode3D
