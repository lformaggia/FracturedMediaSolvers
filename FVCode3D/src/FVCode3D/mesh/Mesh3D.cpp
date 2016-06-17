/*!
 * @file Mesh3D.cpp
 * @brief Classes that implements polyhedrical unstructured meshes (definitions).
 */

#include <FVCode3D/mesh/Mesh3D.hpp>
#include <FVCode3D/geometry/Operations.hpp>
#include <FVCode3D/geometry/BoundingBox.hpp>
#include <FVCode3D/geometry/CoordinateSystem.hpp>

#include <utility>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <array>
#include <queue>

//#define FVCODE3D_DEBUG_CELLS
#define FVCODE3D_CELL_VOLUME_BY_GAUSS

namespace FVCode3D
{

// --------------------   Class Mesh3D::Facet3D   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================
Mesh3D::Facet3D::Facet3D() :
        M_mesh(0), M_area(0.), M_borderId(0), M_zone(0) {}

Mesh3D::Facet3D::Facet3D(const Facet3D & facet) :
        M_mesh(facet.getMesh()), M_vertexIds(facet.getVerticesVector()),
        M_separatedCells(facet.getSeparatedCells()),
        M_representedFractureIds(facet.getRepresentedFractureIds()),
        M_centroid(facet.getCentroid()),
        M_normal(facet.getUnsignedNormal()), M_area(facet.getArea()),
        M_borderId(facet.getBorderId()), M_zone(facet.getZoneCode()) {}

Mesh3D::Facet3D::Facet3D(Mesh3D * const mesh, const std::vector<UInt> & vertices, const UInt zone, const UInt borderId) :
                M_mesh(mesh), M_vertexIds(vertices), M_borderId(borderId), M_zone(zone)
{
    computeCentroidAndNormalAndArea();
}

// ==================================================
// Methods
// ==================================================

void Mesh3D::Facet3D::computeCentroidAndNormalAndArea()
{
    assert(M_vertexIds.size() >= 3);

    if(M_vertexIds.size() > 3)
    {
        Point3D center(0., 0., 0.);

        for(auto id : M_vertexIds)
        {
            center += M_mesh->getNodesVector()[id];
        }
        center /= M_vertexIds.size();

        auto it1 = M_vertexIds.begin();
        auto it2 = it1;
        it2++;

        Real area = 0.;
        M_area = 0.;
        M_centroid = Point3D(0., 0., 0.);
        M_normal = Point3D(0., 0., 0.);

        for( ; it2 != M_vertexIds.end(); ++it1, ++it2)
        {
            const Point3D B = M_mesh->getNodesVector()[*it1];
            const Point3D C = M_mesh->getNodesVector()[*it2];

            area = FVCode3D::triangleArea(center, B, C);

            M_area += area;
            M_centroid += ( center + B + C ) / 3. * area;
            M_normal += FVCode3D::computeNormal(center, B, C) * area;
        }

        const Point3D B = M_mesh->getNodesVector()[ M_vertexIds[M_vertexIds.size()-1] ];
        const Point3D C = M_mesh->getNodesVector()[ M_vertexIds[0] ];

        area = FVCode3D::triangleArea(center, B, C );

        M_area += area;
        M_centroid += ( center + B + C ) / 3. * area;
        M_normal += FVCode3D::computeNormal(center, B, C) * area;

        M_centroid /= M_area;
        M_normal /= M_area;
        M_normal.normalize();
    }
    else
    {
        const Point3D A(M_mesh->getNodesVector()[M_vertexIds[0]]);
        const Point3D B(M_mesh->getNodesVector()[M_vertexIds[1]]);
        const Point3D C(M_mesh->getNodesVector()[M_vertexIds[2]]);
        M_normal = FVCode3D::computeNormal(A,B,C);

        M_area = FVCode3D::triangleArea(A,B,C);

        M_centroid = ( M_mesh->getNodesVector()[M_vertexIds[0]] +
                       M_mesh->getNodesVector()[M_vertexIds[1]] +
                       M_mesh->getNodesVector()[M_vertexIds[2]]
                     ) / 3.;
    }
}

void Mesh3D::Facet3D::computeNormal()
{
    assert(M_vertexIds.size() >= 3);

    if(M_vertexIds.size() > 3)
    {
        Point3D center(0., 0., 0.);

        for(auto id : M_vertexIds)
        {
            center += M_mesh->getNodesVector()[id];
        }
        center /= M_vertexIds.size();

        auto it1 = M_vertexIds.begin();
        auto it2 = it1;
        it2++;
        Real area = 0.;
        M_area = 0.;
        M_normal = Point3D(0., 0., 0.);

        for( ; it2 != M_vertexIds.end(); ++it1, ++it2)
        {
            const Point3D B = M_mesh->getNodesVector()[*it1];
            const Point3D C = M_mesh->getNodesVector()[*it2];

            area = FVCode3D::triangleArea(center, B, C);

            M_area += area;
            M_normal += FVCode3D::computeNormal(center, B, C) * area;
        }

        const Point3D B = M_mesh->getNodesVector()[ M_vertexIds[M_vertexIds.size()-1] ];
        const Point3D C = M_mesh->getNodesVector()[ M_vertexIds[0] ];

        area = FVCode3D::triangleArea( center, B, C );

        M_area += area;
        M_normal += FVCode3D::computeNormal(center, B, C) * area;
        M_normal /= M_area;
        M_normal.normalize();
    }
    else
    {
        const Point3D A(M_mesh->getNodesVector()[M_vertexIds[0]]);
        const Point3D B(M_mesh->getNodesVector()[M_vertexIds[1]]);
        const Point3D C(M_mesh->getNodesVector()[M_vertexIds[2]]);
        M_normal = FVCode3D::computeNormal(A,B,C);
    }

    // --- For planar facet ---
//  UInt i = 3;
//  const UInt nPoints = getNumberOfVertices();
//  const Point3D A(M_mesh->getNodesVector()[M_vertexIds[0]]);
//  const Point3D B(M_mesh->getNodesVector()[M_vertexIds[1]]);
//  Point3D C(M_mesh->getNodesVector()[M_vertexIds[2]]);
//  while( (std::fabs(innerAngleRad(B-A,C-A)) < 1e-6) && (i<nPoints) )
//      C = M_mesh->getNodesVector()[M_vertexIds[i++]];
//
//  M_normal = FVCode3D::computeNormal(A,B,C);
}

void Mesh3D::Facet3D::computeArea()
{
    assert(M_vertexIds.size() >= 3);

    if(M_vertexIds.size() > 3)
    {
        Point3D center(0., 0., 0.);

        for(auto id : M_vertexIds)
        {
            center += M_mesh->getNodesVector()[id];
        }
        center /= M_vertexIds.size();

        auto it1 = M_vertexIds.begin();
        auto it2 = it1;
        it2++;
        M_area = 0.;

        for( ; it2 != M_vertexIds.end(); ++it1, ++it2)
        {
            const Point3D B = M_mesh->getNodesVector()[*it1];
            const Point3D C = M_mesh->getNodesVector()[*it2];

            M_area = FVCode3D::triangleArea(center, B, C);
        }

        const Point3D B = M_mesh->getNodesVector()[ M_vertexIds[M_vertexIds.size()-1] ];
        const Point3D C = M_mesh->getNodesVector()[ M_vertexIds[0] ];

        M_area = FVCode3D::triangleArea( center, B, C );
    }
    else
    {
        const Point3D A(M_mesh->getNodesVector()[M_vertexIds[0]]);
        const Point3D B(M_mesh->getNodesVector()[M_vertexIds[1]]);
        const Point3D C(M_mesh->getNodesVector()[M_vertexIds[2]]);
        M_area = FVCode3D::triangleArea(A,B,C);
    }

    // --- For planar facet, also non-convex ---
//    std::vector<Point3D> points;
//    CoordinateSystem3D cs;
//    Real area(0.);
//    const UInt nPoints = M_vertexIds.size();
//
//    cs.computeCartesianCoordinateSystem(M_normal);
//    points.reserve(nPoints);
//
//    for(UInt i=0; i < nPoints; ++i)
//        points.emplace_back( M_mesh->getNodesVector()[M_vertexIds[i]].convertInLocalCoordinates(cs, Point3D(0.,0.,0.)) );
//
//    for(UInt i=0; i < (nPoints-1) ; ++i)
//        area += points[i].x() * points[i+1].y() - points[i].y() * points[i+1].x();
//
//    area += points[nPoints-1].x() * points[0].y() - points[nPoints-1].y() * points[0].x();
//
//    M_area = std::fabs(area)/2;
}

void Mesh3D::Facet3D::computeCentroid()
{
    assert(M_vertexIds.size() >= 3);

    if(M_vertexIds.size() > 3)
    {
        Point3D center(0., 0., 0.);

        for(auto id : M_vertexIds)
        {
            center += M_mesh->getNodesVector()[id];
        }
        center /= M_vertexIds.size();

        auto it1 = M_vertexIds.begin();
        auto it2 = it1;
        it2++;
        Real area = 0.;
        M_area = 0.;
        M_centroid = Point3D(0., 0., 0.);

        for( ; it2 != M_vertexIds.end(); ++it1, ++it2)
        {
            const Point3D B = M_mesh->getNodesVector()[*it1];
            const Point3D C = M_mesh->getNodesVector()[*it2];

            area = FVCode3D::triangleArea(center, B, C);

            M_area += area;
            M_centroid += ( center + B + C ) / 3. * area;
        }

        const Point3D B = M_mesh->getNodesVector()[ M_vertexIds[M_vertexIds.size()-1] ];
        const Point3D C = M_mesh->getNodesVector()[ M_vertexIds[0] ];

        area = FVCode3D::triangleArea( center, B, C );

        M_area += area;
        M_centroid += ( center + B + C ) / 3. * area;
        M_centroid /= M_area;
    }
    else
    {
        const Point3D A(M_mesh->getNodesVector()[M_vertexIds[0]]);
        const Point3D B(M_mesh->getNodesVector()[M_vertexIds[1]]);
        const Point3D C(M_mesh->getNodesVector()[M_vertexIds[2]]);
        M_centroid = ( M_mesh->getNodesVector()[M_vertexIds[0]] +
                       M_mesh->getNodesVector()[M_vertexIds[1]] +
                       M_mesh->getNodesVector()[M_vertexIds[2]]
                     ) / 3.;
    }

    // --- For planar facet, but not non-convex (use computeArea() for planar facet) ---
//    UInt N = M_vertexIds.size();
//    Real tmpArea, totArea = 0.;
//
//    for (UInt j = 1; j < N-1; ++j)
//    {
//        tmpArea = triangleArea( M_mesh->getNodesVector()[M_vertexIds[0]],
//                                M_mesh->getNodesVector()[M_vertexIds[j]],
//                                M_mesh->getNodesVector()[M_vertexIds[j+1]]);
//        totArea += tmpArea;
//        M_centroid += tmpArea * (   M_mesh->getNodesVector()[M_vertexIds[0]] +
//                                        M_mesh->getNodesVector()[M_vertexIds[j]] +
//                                        M_mesh->getNodesVector()[M_vertexIds[j+1]]
//                                    ) / 3.;
//    }
//
//    M_centroid /= totArea;
}

void Mesh3D::Facet3D::showMe(std::ostream  & out) const
{
    out << "Type = Facet3D :";
    out << "M_mesh = " << M_mesh << std::endl;
    out << "  ->" << std::endl;
    for(UInt i=0; i<M_vertexIds.size(); ++i)
    {
        out << "    M_idP["<<i<<"] = " << M_vertexIds[i] << std::endl;
    }
    out << "  -> M_separatedCells : size = " << M_separatedCells.size();
    out << "    [ ";
    if( !M_separatedCells.empty() )
        for( std::set<UInt>::const_iterator it = M_separatedCells.begin();
                it != M_separatedCells.end(); ++it )
            out << *it << " ";
    out << " ] " << std::endl;
    out << "  -> M_representedFractureIds : [ ";
    if( !M_representedFractureIds.empty() )
        for( std::set<UInt>::const_iterator it = M_representedFractureIds.begin();
                it != M_representedFractureIds.end(); ++it )
            out << *it << " ";
    out << " ] " << std::endl;
}


// --------------------   Class Mesh3D::Cell3D   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================
Mesh3D::Cell3D::Cell3D() :
        M_mesh(0), M_volume(0.), M_zone(0) {}

Mesh3D::Cell3D::Cell3D(const Cell3D & cell) :
        M_mesh(cell.getMesh()), M_vertexIds(cell.getVerticesVector()),
        M_facetIds(cell.getFacetsSet()), M_neighborIds(cell.getNeighborsSet()),
        M_centroid(cell.getCentroid()),
        M_volume(cell.volume()), M_zone(cell.getZoneCode()) {}

Mesh3D::Cell3D::Cell3D( Mesh3D * const mesh, const std::vector<UInt> & facets, const UInt zone ) :
        M_mesh(mesh), M_facetIds(facets.begin(), facets.end()), M_zone(zone)
{
    computeVertexIds();
    computeVolumeAndCentroid();
}

// ==================================================
// Methods
// ==================================================
void Mesh3D::Cell3D::computeVertexIds()
{
    for(std::set<UInt>::const_iterator it=M_facetIds.begin(); it != M_facetIds.end(); ++it)
    {
        const Facet3D & f = M_mesh->getFacetsMap().at(*it);
        const UInt N = f.getNumberOfVertices();
        M_vertexIds.reserve(N);
        for(UInt j=0; j<N; ++j)
            M_vertexIds.push_back(f.getVertexId(j));
    }

    sort(M_vertexIds.begin(), M_vertexIds.end());
    M_vertexIds.erase( unique( M_vertexIds.begin(), M_vertexIds.end() ), M_vertexIds.end() );
}

Point3D Mesh3D::Cell3D::outerNormalToFacet(const UInt & facetId) const throw()
{
    bool found(false);
    Point3D normal(0., 0., 0.), v_centr;

    // control that facetId is a facet of this cell
    for( std::set<UInt>::const_iterator it = M_facetIds.begin(); it != M_facetIds.end(); ++it )
    {
        if( *it == facetId )
            found = true;
    }

    // check if the facet exists in this cell
    if( !found )
    {
        std::stringstream error;
        error << "Error: the facet " << facetId << " is not in the cell.";
        throw std::runtime_error(error.str());
    }

    const Facet3D & facet = M_mesh->getFacetsMap().at( facetId );

    normal = facet.getUnsignedNormal();

    v_centr = M_centroid - facet.getVertex(2);

    if( normal * v_centr > 0 )
        return -normal;

    return normal;
}

bool Mesh3D::Cell3D::hasNeighborsThroughFacet( const UInt & facetId, UInt & idNeighbor) const
{
    bool found(false);

    for( std::set<UInt>::const_iterator it = M_neighborIds.begin(); it != M_neighborIds.end(); ++it )
    {
        if( std::find( M_mesh->getCellsMap().at(*it).getFacetsSet().begin(), M_mesh->getCellsMap().at(*it).getFacetsSet().end(), facetId )
            != M_mesh->getCellsMap().at(*it).getFacetsSet().end() )
        {
            idNeighbor = *it;
            found = 1;
            break;
        }
    }

    return found;
}

void Mesh3D::Cell3D::computeVolumeAndCentroid()
{
    assert(M_vertexIds.size() >= 4);

    if(M_vertexIds.size() > 4)
    {
#ifdef FVCODE3D_CELL_VOLUME_BY_GAUSS
        std::vector<Point3D> nodes;
        std::vector< std::vector<UInt> > facets;
        std::map<UInt, UInt> globalToLocal;

        // Count the number of vertices and facets that
        // define the approximation of the non-planar facets polyhedron
        const UInt nNodes = verticesNumber();
        UInt addedNodes = 0;
        UInt nFacets = 0;
        for(auto id : M_facetIds)
        {
            const UInt nodesFacet = M_mesh->getFacetsMap().at(id).getNumberOfVertices();

            if( nodesFacet > 3 )
            {
                addedNodes++;
                nFacets += nodesFacet;
            }
            else
            {
                nFacets++;
            }
        }

        UInt nodesFacet, i, j, intNodesCount(0), intFacetsCount(0);
        std::set<UInt>::iterator it;

        nodes.resize( nNodes + addedNodes );
        facets.resize( nFacets );

        for(i=0; i < nNodes; ++i)
        {
            nodes[i] = M_mesh->getNodesVector()[M_vertexIds[i]];
            globalToLocal.insert(std::make_pair(M_vertexIds[i], i));
        }

        typedef std::pair< UInt, UInt > edge_Type;
        typedef std::map< edge_Type, std::array<UInt,2> > map_Type;

        std::map< edge_Type, std::array<UInt,2> > edgeToTria;
        std::vector< std::array<edge_Type, 3 > > triaToEdge(nFacets);
        std::pair< map_Type::iterator , bool > itMap;

        for(i=0, it=M_facetIds.begin(); it != M_facetIds.end(); ++i, ++it)
        {
            auto& facet = M_mesh->getFacetsMap().at(*it);
            nodesFacet = facet.getNumberOfVertices();

            if( nodesFacet > 3 )
            {
                Point3D center(0., 0., 0.);

                for(auto id : facet.getVerticesVector())
                {
                    center += M_mesh->getNodesVector()[id];
                }
                center /= nodesFacet;
                nodes[nNodes + intNodesCount] = center;

                // loop over the triangles that settle the current facet
                for(j=0; j < nodesFacet; ++j)
                {
                    // add the j-th triangle
                    facets[intFacetsCount].resize(3);
                    facets[intFacetsCount][0] = nNodes + intNodesCount;
                    facets[intFacetsCount][1] = globalToLocal[ facet.getVerticesVector()[j%nodesFacet] ];
                    facets[intFacetsCount][2] = globalToLocal[ facet.getVerticesVector()[(j+1)%nodesFacet] ];

                    // fill the edgeToTria map and the triaToEdge vector
                    for( UInt e = 0; e < 3; ++e)
                    {
                        itMap = edgeToTria.insert( std::make_pair( std::make_pair(facets[intFacetsCount][e%3],
                                facets[intFacetsCount][(e+1)%3] ), std::array<UInt,2>{ { intFacetsCount, intFacetsCount } } ) );

                        if( !itMap.second )
                            itMap.first->second[1] = intFacetsCount;

                        itMap = edgeToTria.insert( std::make_pair( std::make_pair(facets[intFacetsCount][(e+1)%3],
                                facets[intFacetsCount][e%3] ), std::array<UInt,2>{ { intFacetsCount, intFacetsCount } } ) );

                        if( !itMap.second )
                            itMap.first->second[1] = intFacetsCount;

                        triaToEdge[intFacetsCount][e] = std::make_pair( facets[intFacetsCount][e%3] ,
                                                                        facets[intFacetsCount][(e+1)%3] );
                    }

                    intFacetsCount++;
                }

                intNodesCount++;
            }
            else
            {
                facets[intFacetsCount].resize(3);
                for(j=0; j < 3; ++j)
                    facets[intFacetsCount][j] = globalToLocal[ facet.getVerticesVector()[j] ];

                // fill the edgeToTria map
                for( UInt e = 0; e < 3; ++e)
                {
                    itMap = edgeToTria.insert( std::make_pair( std::make_pair(facets[intFacetsCount][e%3],
                            facets[intFacetsCount][(e+1)%3] ), std::array<UInt,2>{ { intFacetsCount, intFacetsCount } } ) );

                    if( !itMap.second )
                        itMap.first->second[1] = intFacetsCount;

                    itMap = edgeToTria.insert( std::make_pair( std::make_pair(facets[intFacetsCount][(e+1)%3],
                            facets[intFacetsCount][e%3] ), std::array<UInt,2>{ { intFacetsCount, intFacetsCount } } ) );

                    if( !itMap.second )
                        itMap.first->second[1] = intFacetsCount;

                    triaToEdge[intFacetsCount][e] = std::make_pair( facets[intFacetsCount][e%3] ,
                                                                    facets[intFacetsCount][(e+1)%3] );
                }

                intFacetsCount++;
            }
        }

        std::vector<bool> checked(nFacets,false);
        std::vector<bool> orientation(nFacets,true);
        std::queue<UInt> queueTria;
        std::vector<Real> sign(nFacets);

        queueTria.push(0);
        checked[0] = true;
        while(!queueTria.empty())
        {
            const UInt idTria = queueTria.front();
            queueTria.pop();

            for(UInt i = 0; i < 3; ++i)
            {
                const edge_Type edge(facets[idTria][i%3] ,
                                     facets[idTria][(i+1)%3] );
                const auto& trias = *(edgeToTria.find(edge));
                const UInt idOther = ( idTria == trias.second[0] ) ? trias.second[1] : trias.second[0];

                const bool found = std::find( std::begin( triaToEdge[idOther] ) , std::end( triaToEdge[idOther] ) , edge )
                                   != std::end( triaToEdge[idOther] );

                orientation[idOther] = orientation[idTria] ^ found;

                if(!checked[idOther])
                {
                    queueTria.push(idOther);
                    checked[idOther] = true;
                }
            }

            sign[idTria] = ( 2. * static_cast<Real>(orientation[idTria]) - 1. );
        }

#ifdef FVCODE3D_DEBUG_CELLS
        std::cout<<std::endl<<std::endl<<"--- NEW ELEMENT ---"<<std::endl;

        std::cout<<"NODI INIZIALI con centroidi facce"<<std::endl;
        std::cout.precision(16);
        for(i=0; i < nodes.size(); ++i)
        {
            std::cout<<i<<": "<<nodes[i]<<std::endl;
        }

        i=0;
        for(auto& facet : facets)
        {
            std::cout<<"FACCIA "<<i++<<std::endl;
            for(auto& id : facet)
            {
                std::cout<<"\t"<<id;
            }
            std::cout<<std::endl;
        }
//        for(i=0, it=M_facetIds.begin(); it != M_facetIds.end(); ++i, ++it)
//        {
//            auto& facet = M_mesh->getFacetsMap().at(*it);
//            std::cout<<"FACCIA "<<i<<std::endl;
//            for(auto id : facet.getVerticesVector())
//            {
//                std::cout<<"\t"<<globalToLocal[id];
//            }
//            std::cout<<std::endl;
//        }
#endif // FVCODE3D_DEBUG_CELLS

        BoundingBox BB(nodes);
        BB.scaleNodesToUnit(nodes);

#ifdef FVCODE3D_DEBUG_CELLS
        std::cout<<"xMin: "<<BB.xMin()<<std::endl;
        std::cout<<"yMin: "<<BB.yMin()<<std::endl;
        std::cout<<"zMin: "<<BB.zMin()<<std::endl;
        std::cout<<"xMax: "<<BB.xMax()<<std::endl;
        std::cout<<"yMax: "<<BB.yMax()<<std::endl;
        std::cout<<"zMax: "<<BB.zMax()<<std::endl;
        std::cout<<"mx: "<<BB.mx()<<std::endl;
        std::cout<<"my: "<<BB.my()<<std::endl;
        std::cout<<"mz: "<<BB.mz()<<std::endl;
        std::cout<<"qx: "<<BB.qx()<<std::endl;
        std::cout<<"qy: "<<BB.qy()<<std::endl;
        std::cout<<"qz: "<<BB.qz()<<std::endl;

        std::cout<<"NODI SCALATI con centroidi facce"<<std::endl;
        for(i=0; i < nodes.size(); ++i)
        {
            std::cout<<i<<": "<<nodes[i]<<std::endl;
        }
#endif // FVCODE3D_DEBUG_CELLS

        // For each (triangular) facet we build a tetrahedron by adding a
        // fourth point at the center of the cell
        M_volume = 0.;
        M_centroid = Point3D(0., 0., 0.);
        for(auto facet = std::begin( facets ); facet != std::end( facets ); ++facet)
        {
            std::vector<Point3D> p(3);
            p[0] = nodes[(*facet)[0]];
            p[1] = nodes[(*facet)[1]];
            p[2] = nodes[(*facet)[2]];

            const UInt index = std::distance( std::begin( facets ) , facet );
            const Point3D triangleNormal = sign[index] * FVCode3D::computeNormal( p[0], p[1], p[2] );
            const Real triangleArea = FVCode3D::triangleArea( p[0], p[1], p[2] );
            const Point3D triangleCentroid = FVCode3D::triangleCentroid( p[0], p[1], p[2] );
            Point3D localCentroid(0., 0., 0.);

            M_volume += triangleArea * dotProduct( triangleCentroid , triangleNormal );

            localCentroid.x() += ( (p[0] + p[1]).x() * (p[0] + p[1]).x() * triangleNormal.x() +
                                   (p[1] + p[2]).x() * (p[1] + p[2]).x() * triangleNormal.x() +
                                   (p[2] + p[0]).x() * (p[2] + p[0]).x() * triangleNormal.x() );

            localCentroid.y() += ( (p[0] + p[1]).y() * (p[0] + p[1]).y() * triangleNormal.y() +
                                   (p[1] + p[2]).y() * (p[1] + p[2]).y() * triangleNormal.y() +
                                   (p[2] + p[0]).y() * (p[2] + p[0]).y() * triangleNormal.y() );

            localCentroid.z() += ( (p[0] + p[1]).z() * (p[0] + p[1]).z() * triangleNormal.z() +
                                   (p[1] + p[2]).z() * (p[1] + p[2]).z() * triangleNormal.z() +
                                   (p[2] + p[0]).z() * (p[2] + p[0]).z() * triangleNormal.z() );
            localCentroid *= triangleArea / 12.;

            M_centroid += localCentroid;
        }
        M_volume /= 3.;
        M_centroid /= 2. * M_volume;

        M_volume /= (BB.mx() * BB.my() * BB.mz());
        M_volume = std::fabs( M_volume );
        BB.scaleNodesToPhysical(M_centroid);

#ifdef FVCODE3D_DEBUG_CELLS
        std::cout<<"Volume: "<<M_volume<<std::endl;
        std::cout<<"Centroid: "<<M_centroid<<std::endl;
#endif // FVCODE3D_DEBUG_CELLS
//        Real gaussVolume = M_volume;
//        Point3D gaussCentroid = M_centroid;
//        {
#else // FVCODE3D_CELL_VOLUME_BY_GAUSS
        std::vector<Point3D> nodes;
        std::vector< std::vector<UInt> > facets;
        std::map<UInt, UInt> globalToLocal;

        // Count the number of vertices and facets that
        // define the approximation of the non-planar facets polyhedron
        const UInt nNodes = verticesNumber();
        UInt addedNodes = 0;
        UInt nFacets = 0;
        for(auto id : M_facetIds)
        {
            const UInt nodesFacet = M_mesh->getFacetsMap().at(id).getNumberOfVertices();

            if( nodesFacet > 3 )
            {
                addedNodes++;
                nFacets += nodesFacet;
            }
            else
            {
                nFacets++;
            }
        }

        UInt nodesFacet, i, j, intNodesCount(0), intFacetsCount(0);
        std::set<UInt>::iterator it;

        nodes.resize( nNodes + addedNodes );
        facets.resize( nFacets );

        Point3D cellCenter(0., 0., 0.);

        for(i=0; i < nNodes; ++i)
        {
            nodes[i] = M_mesh->getNodesVector()[M_vertexIds[i]];
            globalToLocal.insert(std::make_pair(M_vertexIds[i], i));
            cellCenter += nodes[i];
        }
        cellCenter /= nNodes;

        for(i=0, it=M_facetIds.begin(); it != M_facetIds.end(); ++i, ++it)
        {
            auto& facet = M_mesh->getFacetsMap().at(*it);
            nodesFacet = facet.getNumberOfVertices();

            if( nodesFacet > 3 )
            {
                Point3D center(0., 0., 0.);

                for(auto id : facet.getVerticesVector())
                {
                    center += M_mesh->getNodesVector()[id];
                }
                center /= nodesFacet;
                nodes[nNodes + intNodesCount] = center;

                // loop over the triangles (except the last one) that settle the current facet
                for(j=0; j < nodesFacet - 1; ++j)
                {
                    // add the j-th triangle
                    facets[intFacetsCount].resize(3);
                    facets[intFacetsCount][0] = nNodes + intNodesCount;
                    facets[intFacetsCount][1] = globalToLocal[ facet.getVerticesVector()[j] ];
                    facets[intFacetsCount][2] = globalToLocal[ facet.getVerticesVector()[j+1] ];
                    intFacetsCount++;
                }
                // the last triangle that settles the current facet
                facets[intFacetsCount].resize(3);
                facets[intFacetsCount][0] = nNodes + intNodesCount;
                facets[intFacetsCount][1] = globalToLocal[ facet.getVerticesVector()[nodesFacet - 1] ];
                facets[intFacetsCount][2] = globalToLocal[ facet.getVerticesVector()[0] ];
                intFacetsCount++;

                intNodesCount++;
            }
            else
            {
                facets[intFacetsCount].resize(3);
                for(j=0; j < 3; ++j)
                    facets[intFacetsCount][j] = globalToLocal[ facet.getVerticesVector()[j] ];
                intFacetsCount++;
            }
        }

#ifdef FVCODE3D_DEBUG_CELLS
        std::cout<<std::endl<<std::endl<<"--- NEW ELEMENT ---"<<std::endl;

        std::cout<<"NODI INIZIALI con centroidi facce"<<std::endl;
        std::cout.precision(16);
        for(i=0; i < nodes.size(); ++i)
        {
            std::cout<<i<<": "<<nodes[i]<<std::endl;
        }

        i=0;
        for(auto& facet : facets)
        {
            std::cout<<"FACCIA "<<i++<<std::endl;
            for(auto& id : facet)
            {
                std::cout<<"\t"<<id;
            }
            std::cout<<std::endl;
        }
//        for(i=0, it=M_facetIds.begin(); it != M_facetIds.end(); ++i, ++it)
//        {
//            auto& facet = M_mesh->getFacetsMap().at(*it);
//            std::cout<<"FACCIA "<<i<<std::endl;
//            for(auto id : facet.getVerticesVector())
//            {
//                std::cout<<"\t"<<globalToLocal[id];
//            }
//            std::cout<<std::endl;
//        }
#endif // FVCODE3D_DEBUG_CELLS

        BoundingBox BB(nodes);
        BB.scaleNodesToUnit(nodes);
        BB.scaleNodesToUnit(cellCenter);

#ifdef FVCODE3D_DEBUG_CELLS
        std::cout<<"xMin: "<<BB.xMin()<<std::endl;
        std::cout<<"yMin: "<<BB.yMin()<<std::endl;
        std::cout<<"zMin: "<<BB.zMin()<<std::endl;
        std::cout<<"xMax: "<<BB.xMax()<<std::endl;
        std::cout<<"yMax: "<<BB.yMax()<<std::endl;
        std::cout<<"zMax: "<<BB.zMax()<<std::endl;
        std::cout<<"mx: "<<BB.mx()<<std::endl;
        std::cout<<"my: "<<BB.my()<<std::endl;
        std::cout<<"mz: "<<BB.mz()<<std::endl;
        std::cout<<"qx: "<<BB.qx()<<std::endl;
        std::cout<<"qy: "<<BB.qy()<<std::endl;
        std::cout<<"qz: "<<BB.qz()<<std::endl;

        std::cout<<"NODI SCALATI con centroidi facce"<<std::endl;
        for(i=0; i < nodes.size(); ++i)
        {
            std::cout<<i<<": "<<nodes[i]<<std::endl;
        }
#endif // FVCODE3D_DEBUG_CELLS

        // For each (triangular) facet we build a tetrahedron by adding a
        // fourth point at the center of the cell
        M_volume = 0;
        M_centroid = Point3D(0., 0., 0.);
        for (std::vector< std::vector<UInt> >::const_iterator facetIt = facets.begin();
            facetIt != facets.end(); ++facetIt)
        {
            std::vector<Point3D> tetrahedronNodes(4);
            tetrahedronNodes[0] = nodes[(*facetIt)[0]];
            tetrahedronNodes[1] = nodes[(*facetIt)[1]];
            tetrahedronNodes[2] = nodes[(*facetIt)[2]];
            tetrahedronNodes[3] = cellCenter;
            Real partialVolume = FVCode3D::tetrahedronVolume(tetrahedronNodes);
            M_volume += partialVolume;
            M_centroid += partialVolume * ( tetrahedronNodes[0] + tetrahedronNodes[1] + tetrahedronNodes[2] + tetrahedronNodes[3] ) / 4.;
        }
        M_centroid /= M_volume;

        M_volume /= (BB.mx() * BB.my() * BB.mz());
        BB.scaleNodesToPhysical(M_centroid);

#ifdef FVCODE3D_DEBUG_CELLS
        std::cout<<"Volume: "<<M_volume<<std::endl;
        std::cout<<"Centroid: "<<M_centroid<<std::endl;
#endif // FVCODE3D_DEBUG_CELLS
//        std::cout<<"Volume: "<<M_volume<<std::endl;
//        std::cout<<"GaussVolume: "<<gaussVolume<<std::endl;
//        std::cout<<"Centroid: "<<M_centroid<<std::endl;
//        std::cout<<"GaussCentroid: "<<gaussCentroid<<std::endl;
//
//        std::cout<<"Err Volume: "<<(M_volume - gaussVolume) / M_volume<<std::endl;
//        std::cout<<"Err Centrois: "<<(M_centroid - gaussCentroid).norm() / M_centroid.norm()<<std::endl<<std::endl;
//        M_volume = gaussVolume;
//        M_centroid = gaussCentroid;
//    }
#endif // FVCODE3D_CELL_VOLUME_BY_GAUSS
    }
    else
    {
        std::vector<Point3D> nodes(4);
        nodes[0] = M_mesh->getNodesVector()[M_vertexIds[0]];
        nodes[1] = M_mesh->getNodesVector()[M_vertexIds[1]];
        nodes[2] = M_mesh->getNodesVector()[M_vertexIds[2]];
        nodes[3] = M_mesh->getNodesVector()[M_vertexIds[3]];

        M_volume = FVCode3D::tetrahedronVolume(nodes);
        M_centroid = ( nodes[0] + nodes[1] + nodes[2] + nodes[3] ) / 4.;
    }
}

void Mesh3D::Cell3D::showMe(std::ostream & out) const
{
    out << "Type = Cell3D :";
    out << " M_mesh : " << M_mesh << std::endl;
    out << "  -> M_vertexIds : size = " << M_vertexIds.size() << std::endl;
    for( std::vector<UInt>::const_iterator it = M_vertexIds.begin(); it != M_vertexIds.end(); ++it )
        out << "    VertexId = "<< *it << std::endl;
    out << "  -> M_facetIds : size = " << M_facetIds.size() << std::endl;
    for( std::set<UInt>::const_iterator it = M_facetIds.begin(); it != M_facetIds.end(); ++it )
        out << "    FacetId = "<< *it << std::endl;
}


// --------------------   Class Mesh3D   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================
Mesh3D::Mesh3D():
        M_fn(*this) {}

Mesh3D::Mesh3D(const Mesh3D & mesh):
        M_fn(mesh.getFN()), M_nodes(mesh.getNodesVector()),
        M_facets(mesh.getFacetsMap()), M_cells(mesh.getCellsMap()) {}

void Mesh3D::updateFacetsWithFractures()
{
    for(std::vector<Fracture3D>::iterator it = M_fn.getNetwork().begin(); it != M_fn.getNetwork().end(); ++it)
        for(std::vector<UInt>::iterator jt = it->getFractureFacetsId().begin(); jt != it->getFractureFacetsId().end(); ++jt)
            M_facets[*jt].getRepresentedFractureIds().insert(it->getId());
}

void Mesh3D::updateFacetsWithCells()
{
    for(std::map<UInt,Cell3D>::iterator it = M_cells.begin(); it != M_cells.end(); ++it)
        for(std::set<UInt>::iterator jt = it->second.getFacetsSet().begin(); jt != it->second.getFacetsSet().end(); ++jt)
            M_facets[*jt].getSeparatedCells().insert(it->first);
}

void Mesh3D::updateCellsWithNeighbors()
{
    for(std::map<UInt,Facet3D>::iterator it = M_facets.begin(); it != M_facets.end(); ++it)
    {
        for(std::set<UInt>::iterator jt = it->second.getSeparatedCells().begin(); jt != it->second.getSeparatedCells().end(); ++jt)
        {
            for(std::set<UInt>::reverse_iterator kt = it->second.getSeparatedCells().rbegin(); kt.base() != jt; ++kt)
            {
                M_cells[*jt].getNeighborsSet().insert(*kt);
                M_cells[*kt].getNeighborsSet().insert(*jt);
            }
        }
    }
}

void Mesh3D::buildNodesToFacetMap()
{
    std::map<UInt,Facet3D>::const_iterator itF;
    std::vector<UInt> nodes;

    for(itF = M_facets.begin(); itF != M_facets.end(); ++itF)
    {
        nodes = itF->second.getVerticesVector();
        sort(nodes.begin(), nodes.end());
        M_nodesToFacet.insert( std::pair<std::vector<UInt>, UInt>(nodes, itF->first) );
    }
}

UInt Mesh3D::getFacetFromNodes(std::vector<UInt> & nodes) throw()
{
    bool found = true;
    UInt idFacet = 0;

    if(M_nodesToFacet.empty())
    {
        std::set<UInt> nodesIds(nodes.begin(), nodes.end());
        std::set<UInt> facetsIds;
        std::map<UInt,Facet3D>::const_iterator itF;
        std::set<UInt>::const_iterator it, it2;

        for(itF = M_facets.begin(); itF != M_facets.end(); ++itF)
        {
            found = true;
            idFacet = itF->first;
            facetsIds.clear();
            facetsIds.insert(itF->second.getVerticesVector().begin(), itF->second.getVerticesVector().end());
            for(it = facetsIds.begin(), it2 = nodesIds.begin(); it != facetsIds.end(); ++it, ++it2)
            {
                if( *it != *it2 )
                {
                    found = false;
                    break;
                }
            }
            if(found)
                break;
        }
    }
    else
    {
        std::map<std::vector<UInt>, UInt>::const_iterator itM;
        sort(nodes.begin(), nodes.end());
        itM = M_nodesToFacet.find(nodes);
        if(itM == M_nodesToFacet.end())
            found = false;
        else
            idFacet = itM->second;
    }

    if(!found)
    {
        throw std::runtime_error("Error: facet not found in " + std::string(__FUNCTION__) + "." );
    }

    return idFacet;
}

bool Mesh3D::exportVTU(const std::string & filename) const throw()
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

    UInt nPoints = M_nodes.size();
    UInt nCells = M_cells.size();
    UInt nVal = 0, offsets = 0, faceOffsets = 0;

    for( std::map<UInt,Cell3D>::const_iterator it = M_cells.begin(); it != M_cells.end(); ++it )
        nVal += it->second.verticesNumber() + 1;

    // Header
    filestr << "<?xml version=\"1.0\"?>" << std::endl;
    filestr << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl;
    filestr << "\t<UnstructuredGrid>" << std::endl;
    filestr << "\t\t<Piece NumberOfPoints=\"" << nPoints << "\" NumberOfCells=\"" << nCells << "\">" << std::endl;

    // CellData
    filestr << "\t\t\t<CellData Scalars=\"cellID\">" << std::endl;
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"cellID\" format=\"ascii\">" << std::endl;
    filestr << std::scientific << std::setprecision(0);
    for( std::map<UInt,Cell3D>::const_iterator it = M_cells.begin(); it != M_cells.end(); ++it )
        filestr << it->first << std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;
    filestr << "\t\t\t</CellData>" << std::endl;

    filestr << std::scientific << std::setprecision(16);

    // Points
    filestr << "\t\t\t<Points>" << std::endl;
    filestr << "\t\t\t\t<DataArray type=\"" << "Float64" << "\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
    for( std::vector<Point3D>::const_iterator it = M_nodes.begin(); it != M_nodes.end(); ++it )
        filestr << it->x() << " " << it->y() << " " << it->z() <<std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;
    filestr << "\t\t\t</Points>" << std::endl;

    // Cells
    filestr << "\t\t\t<Cells>" << std::endl;
    //  Connectivity
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
    for( std::map<UInt,Cell3D>::const_iterator it = M_cells.begin(); it != M_cells.end(); ++it )
    {
        for( std::vector<UInt>::const_iterator jt = it->second.getVerticesVector().begin(); jt != it->second.getVerticesVector().end()-1; ++jt )
            filestr << *jt << " ";
        filestr << *(it->second.getVerticesVector().rbegin()) << std::endl;
    }
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    //  Offsets
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << std::endl;
    for( std::map<UInt,Cell3D>::const_iterator it = M_cells.begin(); it != M_cells.end(); ++it )
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
    for( std::map<UInt,Cell3D>::const_iterator it = M_cells.begin(); it != M_cells.end(); ++it )
    {
        filestr << it->second.facetsNumber() << std::endl;
        for( std::set<UInt>::const_iterator jt = it->second.getFacetsSet().begin(); jt != it->second.getFacetsSet().end(); ++jt )
        {
            filestr << M_facets.at(*jt).getNumberOfVertices() << std::endl;
            for( std::vector<UInt>::const_iterator kt = M_facets.at(*jt).getVerticesVector().begin(); kt != M_facets.at(*jt).getVerticesVector().end()-1; ++kt )
                filestr << *kt << " ";
            filestr << *(M_facets.at(*jt).getVerticesVector().rbegin()) << std::endl;
        }
    }
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    //  Faceoffsets
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"faceoffsets\" format=\"ascii\">" << std::endl;
    for( std::map<UInt,Cell3D>::const_iterator it = M_cells.begin(); it != M_cells.end(); ++it )
    {
        faceOffsets += 1 + it->second.facetsNumber();
        for( std::set<UInt>::const_iterator jt = it->second.getFacetsSet().begin(); jt != it->second.getFacetsSet().end(); ++jt )
            faceOffsets += M_facets.at(*jt).getNumberOfVertices();
        filestr << faceOffsets << std::endl;
    }
    filestr << "\t\t\t\t</DataArray>" << std::endl;
    filestr << "\t\t\t</Cells>" << std::endl;

    filestr << "\t\t</Piece>" << std::endl;
    filestr << "\t</UnstructuredGrid>" << std::endl;
    filestr << "\t</VTKFile>" << std::endl;

    filestr.close();

    return true;
}

bool Mesh3D::exportCellsVTU(const std::string & filename, const std::vector<UInt> & idCells) const throw()
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

    UInt nPoints = 0;
    UInt nCells = idCells.size();
    UInt offsets = 0, faceOffsets = 0;

    for( std::vector<UInt>::const_iterator it = idCells.begin(); it != idCells.end(); ++it )
        nPoints += M_cells.at(*it).verticesNumber();

    // Header
    filestr << "<?xml version=\"1.0\"?>" << std::endl;
    filestr << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl;
    filestr << "\t<UnstructuredGrid>" << std::endl;
    filestr << "\t\t<Piece NumberOfPoints=\"" << nPoints << "\" NumberOfCells=\"" << nCells << "\">" << std::endl;

    // CellData
    filestr << "\t\t\t<CellData Scalars=\"cellID\">" << std::endl;
    filestr << "<\t\t\t\tDataArray type=\"Int64\" Name=\"cellID\" format=\"ascii\">" << std::endl;
    filestr << std::scientific << std::setprecision(0);
    for( std::vector<UInt>::const_iterator it = idCells.begin(); it != idCells.end(); ++it )
        filestr << *it << std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;
    filestr << "\t\t\t</CellData>" << std::endl;

    filestr << std::scientific << std::setprecision(16);

    // Pointdata
    filestr << "\t\t\t<Points>" << std::endl;
    filestr << "\t\t\t\t<DataArray type=\"" << "Float64" << "\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
    for( std::vector<UInt>::const_iterator it = idCells.begin(); it != idCells.end(); ++it )
        for( std::vector<UInt>::const_iterator jt = M_cells.at(*it).getVerticesVector().begin(); jt != M_cells.at(*it).getVerticesVector().end(); ++jt)
            filestr << M_nodes[*jt].x() << " " << M_nodes[*jt].y() << " " << M_nodes[*jt].z() <<std::endl;
    filestr << "\t\t\t\t</DataArray>" << std::endl;
    filestr << "\t\t\t</Points>" << std::endl;

    // Cells
    filestr << "\t\t\t<Cells>" << std::endl;
    //  Connectivity
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
    for( std::vector<UInt>::const_iterator it = idCells.begin(); it != idCells.end(); ++it )
    {
        for( std::vector<UInt>::const_iterator jt = M_cells.at(*it).getVerticesVector().begin(); jt != M_cells.at(*it).getVerticesVector().end()-1; ++jt )
            filestr << *jt << " ";
        filestr << *(M_cells.at(*it).getVerticesVector().rbegin()) << std::endl;
    }
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    //  Offsets
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << std::endl;
    for( std::vector<UInt>::const_iterator it = idCells.begin(); it != idCells.end(); ++it )
    {
        offsets += M_cells.at(*it).verticesNumber();
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
    for( std::vector<UInt>::const_iterator it = idCells.begin(); it != idCells.end(); ++it )
    {
        filestr << M_cells.at(*it).facetsNumber() << std::endl;
        for( std::set<UInt>::const_iterator jt = M_cells.at(*it).getFacetsSet().begin(); jt != M_cells.at(*it).getFacetsSet().end(); ++jt )
        {
            filestr << M_facets.at(*jt).getNumberOfVertices() << std::endl;
            for( std::vector<UInt>::const_iterator kt = M_facets.at(*jt).getVerticesVector().begin(); kt != M_facets.at(*jt).getVerticesVector().end()-1; ++kt )
                filestr << *kt << " ";
            filestr << *(M_facets.at(*jt).getVerticesVector().rbegin()) << std::endl;
        }
    }
    filestr << "\t\t\t\t</DataArray>" << std::endl;

    //  Faceoffsets
    filestr << "\t\t\t\t<DataArray type=\"Int64\" Name=\"faceoffsets\" format=\"ascii\">" << std::endl;
    for( std::vector<UInt>::const_iterator it = idCells.begin(); it != idCells.end(); ++it )
    {
        faceOffsets += 1 + M_cells.at(*it).facetsNumber();
        for( std::set<UInt>::const_iterator jt = M_cells.at(*it).getFacetsSet().begin(); jt != M_cells.at(*it).getFacetsSet().end(); ++jt )
            faceOffsets += M_facets.at(*jt).getNumberOfVertices();
        filestr << faceOffsets << std::endl;
    }
    filestr << "\t\t\t\t</DataArray>" << std::endl;
    filestr << "\t\t\t</Cells>" << std::endl;

    filestr << "\t\t</Piece>" << std::endl;
    filestr << "\t</UnstructuredGrid>" << std::endl;
    filestr << "\t</VTKFile>" << std::endl;

    filestr.close();

    return true;
}

bool Mesh3D::exportFractureNetworkVTK(const std::string & filename) const throw()
{
    return M_fn.exportNetworkVTK(filename);
}

bool Mesh3D::exportFractureVTK( const std::string & filename, const UInt & f) const throw()
{
    return M_fn.getFracture(f).exportVTK(filename);
}

void Mesh3D::clear()
{
    M_fn.clear();
    M_nodes.clear();
    M_facets.clear();
    M_cells.clear();
    M_nodesToFacet.clear();
}

// --------------------   Overloading operator< (for Facet3D)  --------------------

bool operator<(const Mesh3D::Facet3D & f1, const Mesh3D::Facet3D & f2)
{
    std::vector<UInt> idF1, idF2;
    UInt sizeF1, sizeF2, minSize;

    idF1 = f1.getVerticesVector();
    std::sort(idF1.begin(), idF1.end());
    idF2 = f2.getVerticesVector();
    std::sort(idF2.begin(), idF2.end());
    sizeF1 = idF1.size();
    sizeF2 = idF2.size();
    minSize = std::min(sizeF1,sizeF2);
    for(UInt i=0; i<minSize; ++i)
    {
        if(idF1[i] < idF2[i])
            return true;
        if(idF1[i] > idF2[i])
            return false;
    }

    return sizeF1<sizeF2;
}

} // namespace FVCode3D
