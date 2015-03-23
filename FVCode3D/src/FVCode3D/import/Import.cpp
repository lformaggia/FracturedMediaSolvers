/*!
 * @file Import.cpp
 * @brief Classes for loading files (definitions).
 */

#include <FVCode3D/mesh/Mesh3D.hpp>
#include <FVCode3D/property/Properties.hpp>
#include <FVCode3D/import/Import.hpp>
#include <FVCode3D/property/Permeability.hpp>

#include <utility>
#include <map>
#include <random>

namespace FVCode3D
{

void Importer::extractBC(const Real theta)
{
    Point3D normal, center, centerFace;
    Real max;
    UInt compMax;

    for(std::map<UInt, Mesh3D::Facet3D>::iterator it = M_mesh.getFacetsMap().begin(); it != M_mesh.getFacetsMap().end(); ++it)
    {
        if (it->second.getSeparatedCells().size() == 1)
        {
            center = M_mesh.getCellsMap().at(*(it->second.getSeparatedCells().begin())).getCentroid();
            centerFace = it->second.getCentroid();

            normal = it->second.getUnsignedNormal();
            center -= centerFace;
            center.normalize();

            if(normal * center > 0.)
                normal = -normal;

            normal = rotateOf(normal, theta);

            max = std::max(std::fabs(normal.x()), std::fabs(normal.y()));
            compMax = std::fabs(normal.x()) > std::fabs(normal.y()) ? 0 : 1;

            max = std::max(max, std::fabs(normal.z()));
            compMax = max > std::fabs(normal.z()) ? compMax : 2;

            if(normal[compMax]<0.)
                it->second.setBorderID(2*compMax+1);
            else
                it->second.setBorderID(2*compMax+2);
        }
    }
}

void Importer::addBCAndFractures(const Real theta)
{
    extractBC(theta);
    addFractures();
}

void ImporterMedit::import(bool fracturesOn) throw()
{
    std::ifstream file;
    file.open(M_filename.c_str(), std::ios_base::in);
    if(!file)
    {
        throw std::runtime_error("Error: file not opened when importing Medit mesh.");
    }

    UInt nNodes, nFacets, nCells;
    UInt zone, maxZone, bcId;
    UInt i, j;
    Real x, y, z;
    std::vector<UInt> tmp, tmpNodes(3), tmpFacets(4);
    std::set<UInt> zones;
    std::string buffer = "";

    std::vector<Point3D> & nodesRef = M_mesh.getNodesVector();
    std::map<UInt, Mesh3D::Facet3D> & facetsRef = M_mesh.getFacetsMap();
    std::map<UInt, Mesh3D::Cell3D> & cellsRef = M_mesh.getCellsMap();

    Properties prop;

    // Read nodes
    const std::string s2findN = "Vertices";
    while(buffer!=s2findN)
        getline(file, buffer);

    file >> nNodes;
    nodesRef.reserve(nNodes);

    for(i=0; i < nNodes; ++i)
    {
        file >> x; file >> y; file >> z;
        file >> buffer;
        nodesRef.emplace_back(x,y,z); // Point3D
    }

    // Read facets
    const std::string s2findT = "Triangles";
    while(buffer!=s2findT)
        getline(file, buffer);

    file >> nFacets;

    std::shared_ptr<PermeabilityScalar> permPtr(new PermeabilityScalar);
    permPtr->setPermeability(1. , 0);

    tmp.resize(3);
    for(i=0; i < nFacets; ++i)
    {
        for(UInt j=0; j < 3; ++j)
        {
            file >> tmp[j];
            tmp[j]--;
        }
        file >> zone;
        bcId = (zone <= 1000 && zone > 0) ? 1 : 0;
        zone = (zone > 1000) && (zone <= 2000) ? zone : 0;
        facetsRef.emplace( std::piecewise_construct, std::forward_as_tuple(i), std::forward_as_tuple(&M_mesh, tmp, (zone)*static_cast<UInt>(fracturesOn), bcId) );
        if(zone > 1 && fracturesOn && zones.find(zone) == zones.end())
        {
            prop.setProperties(1., 1., permPtr);
            M_properties.setZone(zone, prop);
            zones.insert(zone);
        }
    }
    tmp.clear();

    // Read cells
    const std::string s2findE = "Tetrahedra";
    while(buffer!=s2findE)
        getline(file, buffer);

    file >> nCells;

    maxZone = (zones.size() > 0) ? *std::max_element(zones.begin(), zones.end()) : 0;
    prop.setProperties(1., 1., permPtr);
    M_properties.setZone(maxZone+1, prop);

    M_mesh.buildNodesToFacetMap();

    tmp.resize(4);
    for(i=0; i < nCells; ++i)
    {
        // get the nodes that define the cell
        for(UInt j=0; j < 4; ++j)
        {
            file >> tmp[j];
            tmp[j]--;
        }
        file >> buffer;

        // costruisco le facce della cella partendo dai nodi
        for(j=0; j < 4; ++j)
        {
            tmpNodes[0] = tmp[j % 4];
            tmpNodes[1] = tmp[(j+1) % 4];
            tmpNodes[2] = tmp[(j+2) % 4];
            tmpFacets[j] = M_mesh.getFacetFromNodes(tmpNodes);
        }

        cellsRef.emplace( std::piecewise_construct, std::forward_as_tuple(i), std::forward_as_tuple(&M_mesh, tmpFacets, maxZone+1) );
    }
    tmp.clear();

    file.close();
}

void ImporterMedit::addFractures()
{
    FractureNetwork3D FN(M_mesh);
    std::vector<Fracture3D> fracturesVector;
    std::map<UInt, Fracture3D> fracturesMap;
    std::map<UInt, Fracture3D>::iterator itF;
    Point3D normal, center, centerFace;

    for(std::map<UInt, Mesh3D::Facet3D>::iterator it = M_mesh.getFacetsMap().begin(); it != M_mesh.getFacetsMap().end(); ++it)
    {
        if (it->second.getZoneCode() > 0 && it->second.getBorderId()==0)
        {
            itF = fracturesMap.find(it->second.getZoneCode());
            if (itF != fracturesMap.end())
                itF->second.push_back(it->first);
            else
            {
                fracturesMap.emplace(std::piecewise_construct, std::forward_as_tuple(it->second.getZoneCode()), std::forward_as_tuple(M_mesh) );
                fracturesMap.at(it->second.getZoneCode()).push_back(it->first);
                fracturesMap.at(it->second.getZoneCode()).getId() = it->second.getZoneCode();
            }
        }
    }

    fracturesVector.reserve(fracturesMap.size());
    for(itF = fracturesMap.begin();  itF != fracturesMap.end(); ++itF)
        fracturesVector.push_back(itF->second);

    FN.addFractures(fracturesVector);

    M_mesh.addFractureNetwork(FN);
}

void ImporterTetGen::import(bool fracturesOn) throw()
{
    std::ifstream file;
    file.open( (M_filename + ".node").c_str(), std::ios_base::in);
    if(!file)
    {
        throw std::runtime_error("Error: file not opened when importing TetGen mesh (.node).");
    }
    file.close();
    file.open( (M_filename + ".face").c_str(), std::ios_base::in);
    if(!file)
    {
        throw std::runtime_error("Error: file not opened when importing TetGen mesh (.face).");
    }
    file.close();
    file.open( (M_filename + ".ele").c_str(), std::ios_base::in);
    if(!file)
    {
        throw std::runtime_error("Error: file not opened when importing TetGen mesh (.ele).");
    }
    file.close();

    UInt nNodes, nFacets, nCells;
    UInt zone, maxZone, bcId;
    UInt i, j;
    Real x, y, z;
    std::vector<UInt> tmp, tmpNodes(3), tmpFacets(4);
    std::set<UInt> zones;
    std::string buffer = "";

    std::vector<Point3D> & nodesRef = M_mesh.getNodesVector();
    std::map<UInt, Mesh3D::Facet3D> & facetsRef = M_mesh.getFacetsMap();
    std::map<UInt, Mesh3D::Cell3D> & cellsRef = M_mesh.getCellsMap();

    Properties prop;

    // Read nodes
    file.open( (M_filename + ".node").c_str(), std::ios_base::in);

    file >> nNodes;
    nodesRef.reserve(nNodes);
    file >> buffer; // garbage
    file >> buffer; // garbage
    file >> buffer; // garbage

    for(i=0; i < nNodes; ++i)
    {
        file >> buffer; // id, garbage
        file >> x; file >> y; file >> z;
        nodesRef.emplace_back(x,y,z); // Point3D
    }

    file.close();

    // Read facets
    file.open( (M_filename + ".face").c_str(), std::ios_base::in);

    file >> nFacets;
    file >> buffer; // garbage

    std::shared_ptr<PermeabilityScalar> permPtr(new PermeabilityScalar);
    permPtr->setPermeability(1. , 0);

    tmp.resize(3);
    for(i=0; i < nFacets; ++i)
    {
        file >> buffer; // id, garbage
        for(UInt j=0; j < 3; ++j)
        {
            file >> tmp[j];
            tmp[j]--;
        }
        file >> zone;
        bcId = (zone <= 1000 && zone > 0) ? 1 : 0;
        zone = (zone > 1000) && (zone <= 2000) ? zone : 0;
        facetsRef.emplace( std::piecewise_construct, std::forward_as_tuple(i), std::forward_as_tuple(&M_mesh, tmp, (zone)*static_cast<UInt>(fracturesOn), bcId) );
        if(zone > 1 && fracturesOn && zones.find(zone) == zones.end())
        {
            prop.setProperties(1., 1., permPtr);
            M_properties.setZone(zone, prop);
            zones.insert(zone);
        }
    }
    tmp.clear();

    file.close();

    // Read cells
    file.open( (M_filename + ".ele").c_str(), std::ios_base::in);

    file >> nCells;
    file >> buffer; // garbage
    file >> buffer; // garbage

    maxZone = (zones.size() > 0) ? *std::max_element(zones.begin(), zones.end()) : 0;
    prop.setProperties(1., 1., permPtr);
    M_properties.setZone(maxZone+1, prop);

    M_mesh.buildNodesToFacetMap();

    tmp.resize(4);
    for(i=0; i < nCells; ++i)
    {
        file >> buffer; // id, garbage
        // get the nodes that define the cell
        for(UInt j=0; j < 4; ++j)
        {
            file >> tmp[j];
            tmp[j]--;
        }

        // costruisco le facce della cella partendo dai nodi
        for(j=0; j < 4; ++j)
        {
            tmpNodes[0] = tmp[j % 4];
            tmpNodes[1] = tmp[(j+1) % 4];
            tmpNodes[2] = tmp[(j+2) % 4];
            tmpFacets[j] = M_mesh.getFacetFromNodes(tmpNodes);
        }

        cellsRef.emplace( std::piecewise_construct, std::forward_as_tuple(i), std::forward_as_tuple(&M_mesh, tmpFacets, maxZone+1) );
    }
    tmp.clear();

    file.close();
}

void ImporterTPFA::import(bool fracturesOn) throw()
{
    std::ifstream file;
    file.open(M_filename.c_str(), std::ios_base::in);
    if(!file)
    {
        throw std::runtime_error("Error: file not opened when importing TPFA mesh.");
    }

    UInt nNodes, nFacets, nCells, nZones, volumeCorrection;
    UInt nodesFacet, facetsCell, bcId=0;
    Int zone;
    UInt i, j;
    Real x, y, z;
    Real aperture, porosity, permeability;
    std::vector<UInt> tmp;

    std::vector<Point3D> & nodesRef = M_mesh.getNodesVector();
    std::map<UInt, Mesh3D::Facet3D> & facetsRef = M_mesh.getFacetsMap();
    std::map<UInt, Mesh3D::Cell3D> & cellsRef = M_mesh.getCellsMap();

    Properties prop;

    // Read header
    file >> nNodes;
    file >> nFacets;
    file >> nCells;
    file >> nZones;
    file >> volumeCorrection;

    // Read nodes
    nodesRef.reserve(nNodes);

    for(i=0; i < nNodes; ++i)
    {
        file >> x; file >> y; file >> z;
        nodesRef.emplace_back(x,y,z); // Point3D
    }

    // Read facets
    for(i=0; i < nFacets; ++i)
    {
        file >> nodesFacet;
        tmp.resize(nodesFacet);

        for(j=0; j < nodesFacet; ++j)
            file >> tmp[j];
        file >> zone;
        facetsRef.emplace( std::piecewise_construct, std::forward_as_tuple(i), std::forward_as_tuple(&M_mesh, tmp, (zone+1)*static_cast<UInt>(fracturesOn), bcId) );
    }

    // Read cells
    for(i=0; i < nCells; ++i)
    {
        file >> facetsCell;
        tmp.resize(facetsCell);

        for(j=0; j < facetsCell; ++j)
            file >> tmp[j];
        file >> zone;
        cellsRef.emplace( std::piecewise_construct, std::forward_as_tuple(i), std::forward_as_tuple(&M_mesh, tmp, zone+1) );
    }
    tmp.clear();

    std::shared_ptr<PermeabilityScalar> permPtr(new PermeabilityScalar);

    // Read properties
    for(i=0; i < nZones; ++i)
    {
        file >> zone;
        file >> aperture;
        file >> porosity;
        file >> permeability; // dummy parameter
        file >> permeability;

        permPtr->setPermeability(permeability, 0);

        prop.setProperties(aperture, porosity, permPtr);
        M_properties.setZone(zone+1, prop);
    }

    file.close();
}

void ImporterTPFA::addFractures()
{
    FractureNetwork3D FN(M_mesh);
    std::vector<Fracture3D> fracturesVector;
    std::map<UInt, Fracture3D> fracturesMap;
    std::map<UInt, Fracture3D>::iterator itF;
    Point3D normal, center, centerFace;

    for(std::map<UInt, Mesh3D::Facet3D>::iterator it = M_mesh.getFacetsMap().begin(); it != M_mesh.getFacetsMap().end(); ++it)
    {
        if (it->second.getZoneCode() > 0 && it->second.getBorderId()==0)
        {
            itF = fracturesMap.find(it->second.getZoneCode());
            if (itF != fracturesMap.end())
                itF->second.push_back(it->first);
            else
            {
                fracturesMap.emplace(std::piecewise_construct, std::forward_as_tuple(it->second.getZoneCode()), std::forward_as_tuple(M_mesh) );
                fracturesMap.at(it->second.getZoneCode()).push_back(it->first);
                fracturesMap.at(it->second.getZoneCode()).getId() = it->second.getZoneCode();
            }
        }
    }

    fracturesVector.reserve(fracturesMap.size());
    for(itF = fracturesMap.begin();  itF != fracturesMap.end(); ++itF)
        fracturesVector.push_back(itF->second);

    FN.addFractures(fracturesVector);

    M_mesh.addFractureNetwork(FN);
}

void ImporterForSolver::import(bool fracturesOn) throw()
{
    std::ifstream file;
    file.open(M_filename.c_str(), std::ios_base::in);
    if(!file)
    {
        throw std::runtime_error("Error: file not opened when importing .fvg mesh.");
    }

    UInt nNodes, nFacets, nCells, nFractures;
    UInt nodesFacet, facetsCell, facetsFracture, nSepCells, bcId, isFrac;
    UInt zone, facetId;
    UInt i, j;
    Real x, y, z;
    Real aperture, porosity, permeability;
    std::vector<UInt> tmp;
    std::string buffer = "";
    const UInt permSize = 6;

    std::vector<Point3D> & nodesRef = M_mesh.getNodesVector();
    std::map<UInt, Mesh3D::Facet3D> & facetsRef = M_mesh.getFacetsMap();
    std::map<UInt, Mesh3D::Cell3D> & cellsRef = M_mesh.getCellsMap();
    FractureNetwork3D & FN = M_mesh.getFN();
    Properties prop;

    // Read nodes
    const std::string s2findN = "POINTS";
    while(buffer!=s2findN)
        file >> buffer;

    file >> nNodes;
    nodesRef.reserve(nNodes);

    for(i=0; i < nNodes; ++i)
    {
        file >> x; file >> y; file >> z;
        nodesRef.emplace_back(x,y,z); // Point3D
    }

    // Read facets with properties
    const std::string s2findF = "FACETS";
    while(buffer!=s2findF)
        file >> buffer;

    file >> nFacets;

    PermPtr_Type permPtr;

    zone = 1;
    for(i=0; i < nFacets; ++i)
    {
        file >> nodesFacet;
        tmp.resize(nodesFacet);

        for(UInt j=0; j < nodesFacet; ++j)
            file >> tmp[j];
        file >> nSepCells;
        for(UInt j=0; j<nSepCells; ++j)
            file >> buffer;
        file >> bcId;
        file >> isFrac;
        if (bcId && isFrac)
        {
            isFrac = 0;
            for(UInt j=0; j<permSize + 2; ++j)
                file >> buffer;
        }

        if(isFrac)
        {
            file >> aperture;
            file >> porosity;
            switch(permSize)
            {
                case 1:
                    permPtr.reset( new PermeabilityScalar);
                    break;
                case 3:
                    permPtr.reset( new PermeabilityDiagonal);
                    break;
                case 6:
                    permPtr.reset( new PermeabilitySymTensor);
                    break;
                case 9:
                    permPtr.reset( new PermeabilityFullTensor);
                    break;
            }

            for(UInt j=0; j<permSize; ++j)
            {
                file >> permeability;
                permPtr->operator[](j) = permeability;
            }
        }
        facetsRef.emplace( std::piecewise_construct, std::forward_as_tuple(i),
                           std::forward_as_tuple(&M_mesh, tmp,
                           isFrac*static_cast<UInt>(fracturesOn)*zone, bcId) );
        if (isFrac*static_cast<UInt>(fracturesOn))
        {
            prop.setProperties(aperture, porosity, permPtr);
            M_properties.setZone(zone, prop);
            zone++;
        }
    }

    // Read cells with properties
    const std::string s2findC = "CELLS";
    while(buffer!=s2findC)
        file >> buffer;

    file >> nCells;

    aperture = 1.;
    for(i=0; i < nCells; ++i)
    {
        file >> facetsCell;
        tmp.resize(facetsCell);

        for(UInt j=0; j < facetsCell; ++j)
            file >> tmp[j];

        file >> porosity;
        switch(permSize)
        {
            case 1:
                permPtr.reset( new PermeabilityScalar);
                break;
            case 3:
                permPtr.reset( new PermeabilityDiagonal);
                break;
            case 6:
                permPtr.reset( new PermeabilitySymTensor);
                break;
            case 9:
                permPtr.reset( new PermeabilityFullTensor);
                break;
        }

        for(UInt j=0; j<permSize; ++j)
        {
            file >> permeability;
            permPtr->operator[](j) = permeability;
        }

        cellsRef.emplace( std::piecewise_construct, std::forward_as_tuple(i), std::forward_as_tuple(&M_mesh, tmp, zone) );
        prop.setProperties(aperture, porosity, permPtr);
        M_properties.setZone(zone, prop);
        zone++;
    }
    tmp.clear();

    // Read fracture network
    const std::string s2findFN= "FRACTURE_NETWORK";
    while(buffer!=s2findFN)
        file >> buffer;

    file >> nFractures;

    FN.getNetwork().reserve(nFractures);

    const UInt fracOn = nFractures*static_cast<UInt>(fracturesOn);
    for(i=0; i < fracOn; ++i)
    {
        file >> facetsFracture;
        tmp.reserve(facetsFracture);
        for(j=0; j < facetsFracture; ++j)
        {
            file >> facetId;
            if( ! facetsRef[facetId].isBorderFacet() )
                tmp.push_back(facetId);
        }
        if(tmp.size()>0)
            FN.getNetwork().emplace_back(M_mesh, tmp, i); // Fracture3D
        tmp.clear();
    }

    file.close();
}

} // namespace FVCode3D
