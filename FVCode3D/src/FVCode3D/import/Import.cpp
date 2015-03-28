/*!
 * @file Import.cpp
 * @brief Classes for loading files (definitions).
 */

#include <FVCode3D/mesh/Mesh3D.hpp>
#include <FVCode3D/property/Properties.hpp>
#include <FVCode3D/import/Import.hpp>
#include <FVCode3D/property/Permeability.hpp>
#include <FVCode3D/utility/StringManipolator.hpp>

#include <utility>
#include <map>
#include <random>
#include <locale>

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

void Importer::addFractures()
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

void ImporterOpenFOAM::import(bool fracturesOn) throw()
{
    std::ifstream file;

    file.open( (M_filename + "points").c_str(), std::ios_base::in);
    if(!file)
    {
        throw std::runtime_error("Error: file not opened when importing points.");
    }
    file.close();
    file.open( (M_filename + "faces").c_str(), std::ios_base::in);
    if(!file)
    {
        throw std::runtime_error("Error: file not opened when importing faces.");
    }
    file.close();
    file.open( (M_filename + "owner").c_str(), std::ios_base::in);
    if(!file)
    {
        throw std::runtime_error("Error: file not opened when importing owner.");
    }
    file.close();
    file.open( (M_filename + "neighbour").c_str(), std::ios_base::in);
    if(!file)
    {
        throw std::runtime_error("Error: file not opened when importing neighbour.");
    }
    file.close();
    file.open( (M_filename + "boundary").c_str(), std::ios_base::in);
    if(!file)
    {
        throw std::runtime_error("Error: file not opened when importing boundary.");
    }
    file.close();

    UInt nNodes, nFacets, nOwner, nNeigh, nPatch;
    UInt nodesFacet;
    UInt i, j;
    Real x,y,z;
    UInt cellId, patchId, maxPatchId = 0;
    std::vector<UInt> tmp;

    enum class faceType : SUInt
    {
        None, Internal, Boundary, Fracture
    };

    std::map<Point3D, UInt> vertexMap;
    std::map<Mesh3D::Facet3D, UInt> facetMap;
    std::pair< std::map<Point3D,UInt>::iterator, bool > itMap;
    std::pair< std::map<Mesh3D::Facet3D,UInt>::iterator, bool > itMapF;
    std::vector<UInt> fullToLocal;
    std::vector<UInt> fullToLocalFaces;
    std::vector<faceType> typeFaces;
    UInt vertexId = 0;
    UInt facetId = 0;

    std::vector<Point3D> & nodesRef = M_mesh.getNodesVector();
    std::map<UInt, Mesh3D::Facet3D> & facetsRef = M_mesh.getFacetsMap();
    std::map<UInt, Mesh3D::Cell3D> & cellsRef = M_mesh.getCellsMap();
    Properties prop;
    std::shared_ptr<PermeabilityScalar> permPtr(new PermeabilityScalar);
    permPtr->setPermeability(1. , 0);

    std::string buffer = "";
    std::vector< std::string > tokens;
    std::vector< std::string > subTokens;
    const std::string initString =
        "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //";

    // Read nodes
    file.open( (M_filename + "points").c_str(), std::ios_base::in);

    while(buffer != initString)
    {
        std::getline(file, buffer);
    }

    file.imbue(std::locale(file.getloc(), new MyDelimiters));

    buffer.clear();
    while(buffer == "")
    {
        file >> buffer;
    }

    nNodes = lexical_cast<UInt>(buffer);
    nodesRef.reserve(nNodes);
    fullToLocal.resize(nNodes);

    for(i=0; i < nNodes; ++i)
    {
        file >> x; file >> y; file >> z;

        Point3D V(x, y, z);
        itMap = vertexMap.insert(std::pair<FVCode3D::Point3D,UInt>(V, 0));

        if (itMap.second)
        {
            itMap.first->second = vertexId;
            nodesRef.emplace_back(x, y, z); // Point3D
            ++vertexId;
        }
        fullToLocal[i] = itMap.first->second;
    }
    nNodes = vertexId;

    file.imbue(std::locale::classic());
    file.close();

    // Read faces
    file.open( (M_filename + "faces").c_str(), std::ios_base::in);

    while(buffer != initString)
    {
        std::getline(file, buffer);
    }

    file.imbue(std::locale(file.getloc(), new MyDelimiters));

    buffer.clear();
    while(buffer == "")
    {
        file >> buffer;
    }

    nFacets = lexical_cast<UInt>(buffer);
    fullToLocalFaces.resize(nFacets);

    for(i=0; i < nFacets; ++i)
    {
        file >> nodesFacet;
        tmp.resize(nodesFacet);

        for(j=0; j < nodesFacet; ++j)
        {
            file >> tmp[j];
            tmp[j] = fullToLocal[tmp[j]];
        }

        Mesh3D::Facet3D facet(&M_mesh, tmp, 0, 0);
        itMapF = facetMap.insert(std::pair<Mesh3D::Facet3D, UInt>(facet,0));

        if(itMapF.second)
        {
            itMapF.first->second = facetId;
            facetsRef.emplace( std::piecewise_construct, std::forward_as_tuple(facetId), std::forward_as_tuple(&M_mesh, tmp, 0, 0) );
            ++facetId;
        }
        fullToLocalFaces[i] = itMapF.first->second;
    }
    nFacets = facetId;

    file.imbue(std::locale::classic());
    file.close();

    // Read owner
    file.open( (M_filename + "owner").c_str(), std::ios_base::in);

    while(buffer != initString)
    {
        std::getline(file, buffer);
    }

    file.imbue(std::locale(file.getloc(), new MyDelimiters));

    buffer = "";
    while(buffer == "")
    {
        file >> buffer;
    }

    nOwner = lexical_cast<UInt>(buffer);
    typeFaces.resize(nFacets, faceType::None);

    for(i=0; i < nOwner; ++i)
    {
        file >> cellId;
        if( facetsRef[fullToLocalFaces[i]].getSeparatedCells().size() > 0 )
        {
            typeFaces[fullToLocalFaces[i]] = faceType::Fracture;
        }
        facetsRef[fullToLocalFaces[i]].getSeparatedCells().insert(cellId);

        cellsRef.emplace( std::piecewise_construct, std::forward_as_tuple(cellId), std::forward_as_tuple() );
        cellsRef[cellId].setMesh(&M_mesh);
        cellsRef[cellId].getFacetsSet().insert(fullToLocalFaces[i]);
    }
    file.imbue(std::locale::classic());
    file.close();

// Probably this is useless. I don't know how the double faces are managed in openFoam.
// Enable this piece of code if you can't load the fractures, and try...
//    for(auto facet : facetsRef)
//    {
//        if( facet.second.getSeparatedCells().size() == 2 )
//        {
//            typeFaces[facet.first] = faceType::Fracture;
//        }
//    }

    // Read neighbour
    file.open( (M_filename + "neighbour").c_str(), std::ios_base::in);

    while(buffer != initString)
    {
        std::getline(file, buffer);
    }

    file.imbue(std::locale(file.getloc(), new MyDelimiters));

    buffer = "";
    while(buffer == "")
    {
        file >> buffer;
    }

    nNeigh = lexical_cast<UInt>(buffer);

    for(i=0; i < nNeigh; ++i)
    {
        file >> cellId;
        facetsRef[fullToLocalFaces[i]].getSeparatedCells().insert(cellId);

        cellsRef.emplace( std::piecewise_construct, std::forward_as_tuple(cellId), std::forward_as_tuple() );
        cellsRef[cellId].setMesh(&M_mesh);
        cellsRef[cellId].getFacetsSet().insert(fullToLocalFaces[i]);
    }
    file.imbue(std::locale::classic());
    file.close();

    for(auto facet : facetsRef)
    {
        if( typeFaces[facet.first] == faceType::None)
        {
            if( facet.second.getSeparatedCells().size() == 1 )
            {
                typeFaces[facet.first] = faceType::Boundary;
            }
            else
            {
                typeFaces[facet.first] = faceType::Internal;
            }
        }
    }

    // Read boundary
    file.open( (M_filename + "boundary").c_str(), std::ios_base::in);

    while(buffer != initString)
    {
        std::getline(file, buffer);
    }

    file.imbue(std::locale(file.getloc(), new MyDelimiters));

    buffer = "";
    while(buffer == "")
    {
        file >> buffer;
    }

    nPatch = lexical_cast<UInt>(buffer);

    for(i=0; i < nPatch; ++i)
    {
        std::size_t found;
        UInt startId, patchLength;

        do
        {
            std::getline(file, buffer);
            found = buffer.find("patch");
        }
        while(found==std::string::npos);

        patchId = lexical_cast<UInt>(buffer.substr(found+5));
        maxPatchId = (patchId > maxPatchId) ? patchId : maxPatchId;

        for(j=0; j<6; ++j)
        {
            file >> buffer; // garbage
        }

        file >> patchLength;
        file >> buffer; // garbage
        file >> startId;

        for(j=startId; j < (startId+patchLength); ++j)
        {
            if(typeFaces[fullToLocalFaces[j]]==faceType::Boundary)
            {
                facetsRef[fullToLocalFaces[j]].setBorderID(patchId);
            }
            else if (fracturesOn)
            {
                facetsRef[fullToLocalFaces[j]].setZoneCode( patchId );
                prop.setProperties(1., 1., permPtr);
                M_properties.setZone(patchId, prop);
            }
        }
    }

    file.imbue(std::locale::classic());
    file.close();

    // build cells
    ++maxPatchId;
    for(auto & cell : cellsRef)
    {
        cell.second.setZoneCode(maxPatchId);
        cell.second.computeVertexIds();
        cell.second.computeVolumeAndCentroid();
    }
    prop.setProperties(1., 1., permPtr);
    M_properties.setZone(maxPatchId, prop);

    // clear the separated cells from the facets
    for(auto & facet : facetsRef)
    {
        facet.second.getSeparatedCells().clear();
    }
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
