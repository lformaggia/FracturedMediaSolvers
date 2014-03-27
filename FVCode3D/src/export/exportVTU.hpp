/*!
 *  @file exportVTU.hpp
 *  @brief Classes for saving files in the VTU format.
 */

#ifndef EXPORTVTU_HPP_
#define EXPORTVTU_HPP_

#include "export/export.hpp"

/*!
 * @class ExporterVTU
 * This class allows to export mesh, fractures, solution, and properties in VTU format.
 */
class ExporterVTU : public Exporter
{
public:

    //! Constructor
    ExporterVTU(){}

    //! Export the mesh (only cells)
    /*!
     * @param mesh reference of a Geometry::Mesh3D
     * @param filename name of the file
     */
    virtual void exportMesh(const Mesh3D & mesh, const std::string filename);

    //! Export the fracture facets
    /*!
     * @param mesh reference of a Geometry::Mesh3D
     * @param filename name of the file
     */
    virtual void exportFractures(const Mesh3D & mesh, const std::string filename);

    //! Export the mesh, cells and fracture facets, in a single file
    /*!
     * @param mesh reference of a Geometry::Mesh3D
     * @param filename name of the file
     */
    virtual void exportMeshWithFractures(const Mesh3D & mesh, const std::string filename);

    //! Export the wireframe
    /*!
     * @param mesh reference of a Geometry::Mesh3D
     * @param filename name of the file
     */
    virtual void exportWireframe(const Mesh3D & mesh, const std::string filename);
    
    //! Export the fracture junctures
    /*!
     * @param mesh reference of a Geometry::Rigid_Mesh
     * @param filename name of the file
     */
    virtual void exportFractureJunctures(const Rigid_Mesh & mesh, const std::string filename);

    //! Export the solution on cells and fracture facets in a single file
    /*!
     * @param mesh reference of a Geometry::Rigid_Mesh
     * @param filename name of the file
     * @param sol Eigen vector that contain the solution (cells + fracture facets)
     */
    template <typename VectorType>
    void exportSolution(const Rigid_Mesh & mesh, const std::string filename, const VectorType & sol);

    //! Export the solution on fracture facets
    /*!
     * @param mesh reference of a Geometry::Rigid_Mesh
     * @param filename name of the file
     * @param sol Eigen vector that contain the solution (cells + fracture facets)
     */
    virtual void exportSolutionOnFractures(const Rigid_Mesh & mesh, const std::string filename, const Eigen::VectorXd & sol);

    //! Export the a specific property on cells and fracture facets
    /*!
     * @param mesh reference of a Geometry::Mesh3D
     * @param properties reference to a Geometry::PropertiesMap
     * @param filename name of the file
     * @param propertiesType flag used to select which properties to export
     * @param property pointer to a generic property
     */
    virtual void exportWithProperties(const Mesh3D & mesh, const PropertiesMap & properties, const std::string filename, const Flag16bit propertiesType, const std::vector<Real> * property = NULL);

    //! Export all properties defined on cells and fracture facets
    /*!
     * @param mesh reference of a Geometry::Mesh3D
     * @param properties reference to a Geometry::PropertiesMap
     * @param filename name of the file
     */
    virtual void exportWithProperties(const Mesh3D & mesh, const PropertiesMap & properties, const std::string filename);

    //! Export the solution on fracture facets
    /*!
     * @param mesh reference of a Geometry::Rigid_Mesh
     * @param filename name of the file
     * @param sol Eigen vector that contain the solution (cells + fracture facets)
     */
    virtual void exportWithProperties(const Rigid_Mesh & mesh, const std::string filename, const Flag16bit propertiesType, const std::vector<Real> * property = NULL ) ;

    //! Destructor
    virtual ~ExporterVTU() {};
};

template <typename VectorType>
void ExporterVTU::exportSolution(const Rigid_Mesh & mesh, const std::string filename, const VectorType & sol)
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
    for( UInt i = 0; i < nTotal; ++i )
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

#endif /* EXPORTVTU_HPP_ */
