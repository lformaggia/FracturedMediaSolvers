/*!
 * @file converter.cpp
 * @brief Methods to convert format files (definitions).
 */

#include <FVCode3D/mesh/Mesh3D.hpp>
#include <FVCode3D/property/Properties.hpp>
#include <FVCode3D/utility/Converter.hpp>
#include <FVCode3D/property/Permeability.hpp>

namespace FVCode3D
{

void saveAsSolverFormat(const std::string filename, Mesh3D & mesh, PropertiesMap & properties) throw()
{
    std::fstream file;

    file.open (filename.c_str(), std::ios_base::out);

    if (file.is_open())
    {
        std::cout << std::endl << " File: " << filename << ", successfully opened";
    }
    else
    {
        throw std::runtime_error("Error: file " + filename + " not opened.");
        return;
    }

    std::cout << std::endl << " Save as solver format... " << std::endl;

    UInt nNodes, nFacets, nCells, nFractures;
    UInt nodesFacet, facetsCell, facetsFracture;
    UInt i, j;

    std::vector<Point3D> & nodesRef = mesh.getNodesVector();
    std::map<UInt, Mesh3D::Facet3D> & facetsRef = mesh.getFacetsMap();
    std::map<UInt, Mesh3D::Cell3D> & cellsRef = mesh.getCellsMap();
    FractureNetwork3D & FN = mesh.getFN();

    nNodes = nodesRef.size();
    nFacets = facetsRef.size();
    nCells = cellsRef.size();
    nFractures = FN.size();

    file << std::scientific << std::setprecision(0);

    // Save header
    file << "# Mesh3D SolverDataFile" << std::endl << std::endl;

    // Save nodes
    file << "POINTS" << std::endl;
    file << nNodes << std::endl;
    file << std::scientific << std::setprecision(10);
    for(i=0; i < nNodes; ++i)
    {
        file << nodesRef[i].x() << " ";
        file << nodesRef[i].y() << " ";
        file << nodesRef[i].z() << std::endl;
    }
    file << std::endl;
    file << std::scientific << std::setprecision(0);

    // Save facets
    file << "FACETS" << std::endl;
    file << nFacets << std::endl;
    for(i=0; i < nFacets; ++i)
    {
        file << std::scientific << std::setprecision(0);
        nodesFacet = facetsRef[i].getNumberOfVertices();
        file << nodesFacet << " ";

        for(j=0; j < nodesFacet; ++j)
            file << facetsRef[i].getVertexId(j) << " ";

        file << facetsRef[i].getSeparatedCells().size() << " ";
        for(std::set<UInt>::const_iterator it = facetsRef[i].getSeparatedCells().begin(); it != facetsRef[i].getSeparatedCells().end(); ++it)
            file << *it << " ";

        file << facetsRef[i].getBorderId() << " ";
        file << facetsRef[i].isFracture() << " ";

        file.unsetf(std::ios_base::floatfield);
        file << /*std::scientific <<*/ std::setprecision(10);
        if(facetsRef[i].isFracture())
        {
            file << properties.getProperties(facetsRef[i].getZoneCode()).M_aperture << " ";
            file << properties.getProperties(facetsRef[i].getZoneCode()).M_porosity;
            for(UInt row = 0; row < 3; ++row)
            {
                for(UInt col = 0; col < 3; ++col)
                {
                    if(col >= row)
                    {
                        file << " "
                             << properties.getProperties(facetsRef[i].getZoneCode()).M_permeability->operator()(row,col);
                    }
                }
            }
        }
        file << std::endl;
    }
    file << std::endl;
    file << std::scientific << std::setprecision(0);

    // Save cells
    file << "CELLS" << std::endl;
    file << nCells << std::endl;
    for(i=0; i < nCells; ++i)
    {
        file << std::scientific << std::setprecision(0);
        facetsCell = cellsRef[i].facetsNumber();
        file << facetsCell << " ";

        for(std::set<UInt>::const_iterator it = cellsRef[i].getFacetsSet().begin(); it != cellsRef[i].getFacetsSet().end(); ++it)
            file << *it << " ";

        file.unsetf(std::ios_base::floatfield);
        file << /*std::scientific <<*/ std::setprecision(10);
        file << properties.getProperties(cellsRef[i].getZoneCode()).M_porosity;
        for(UInt row = 0; row < 3; ++row)
        {
            for(UInt col = 0; col < 3; ++col)
            {
                if(col >= row)
                {
                    file << " "
                         << properties.getProperties(cellsRef[i].getZoneCode()).M_permeability->operator()(row,col);
                }
            }
        }
        file << std::endl;
    }
    file << std::endl;
    file << std::scientific << std::setprecision(0);

    // Save fracture network
    file << "FRACTURE_NETWORK" << std::endl;
    file << nFractures << std::endl;
    for(i=0; i < nFractures; ++i)
    {
        facetsFracture = FN.getFracture(i).getNumberOfFractureFacets();
        file << facetsFracture << " ";

        for(j=0; j < facetsFracture-1; ++j)
            file << FN.getFracture(i).getFractureFacetsId()[j] << " ";
        file << FN.getFracture(i).getFractureFacetsId()[facetsFracture-1] << std::endl;
    }

    file.close();
}

void saveAsMeditFormat(const std::string filename, Mesh3D & mesh) throw()
{
    std::fstream file;

    file.open (filename.c_str(), std::ios_base::out);

    if (file.is_open())
    {
        std::cout << std::endl << " File: " << filename << ", successfully opened";
    }
    else
    {
        throw std::runtime_error("Error: file " + filename + " not opened.");
        return;
    }

    std::cout << std::endl << " Save as Medit format... " << std::endl;

    UInt nNodes, nFacets, nCells;
    UInt i, zone;

    std::vector<Point3D> & nodesRef = mesh.getNodesVector();
    std::map<UInt, Mesh3D::Facet3D> & facetsRef = mesh.getFacetsMap();
    std::map<UInt, Mesh3D::Cell3D> & cellsRef = mesh.getCellsMap();

    nNodes = nodesRef.size();
    nFacets = facetsRef.size();
    nCells = cellsRef.size();

    file << "MeshVersionFormatted 1" << std::endl;
    file << std::endl;
    file << "Dimension" << std::endl << "3" << std::endl;
    file << std::endl;
    file << "Vertices" << std::endl;
    file << nNodes << std::endl;

    for(i=0; i<nNodes; ++i)
    {
        file << nodesRef[i].x() << "  ";
        file << nodesRef[i].y() << "  ";
        file << nodesRef[i].z() << "  ";
        file << "0" << std::endl;
    }

    file << std::endl;
    file << "Triangles" << std::endl;
    file << nFacets << std::endl;
    for(i=0; i<nFacets; ++i)
    {
        file << facetsRef.at(i).getVerticesVector()[0] + 1 << "  ";
        file << facetsRef.at(i).getVerticesVector()[1] + 1 << "  ";
        file << facetsRef.at(i).getVerticesVector()[2] + 1 << "  ";
        zone = facetsRef.at(i).isBorderFacet() ? 1 : 0 ;
        zone = facetsRef.at(i).isFracture() ? 1001 : zone ;
        file << zone << std::endl;
    }

    file << std::endl;
    file << "Tetrahedra" << std::endl;
    file << nCells << std::endl;
    for(i=0; i<nCells; ++i)
    {
        file << cellsRef.at(i).getVerticesVector()[0] + 1 << "  ";
        file << cellsRef.at(i).getVerticesVector()[1] + 1 << "  ";
        file << cellsRef.at(i).getVerticesVector()[2] + 1 << "  ";
        file << cellsRef.at(i).getVerticesVector()[3] + 1 << "  ";
        file << "0" << std::endl;
    }

    file.close();
}

} // namespace FVCode3D
