/*!
 *	@file QuadratureRules.cpp
 *	@brief Classes for quadrature rules (definitions).
 */

#include <FVCode3D/mesh/RigidMesh.hpp>
#include <FVCode3D/quadrature/QuadratureRules.hpp>

namespace FVCode3D
{

std::unique_ptr<QuadratureRule> CentroidQuadrature::clone() const
{
	return std::unique_ptr<QuadratureRule>(new CentroidQuadrature(*this));
}

Real CentroidQuadrature::apply(const Cell & cell, const std::function<Real(Point3D)> & integrand) const
{
	return integrand(cell.getCentroid()) * cell.getVolume();
}

Real CentroidQuadrature::apply(const Facet & facet, const std::function<Real(Point3D)> & integrand) const
{
	return integrand(facet.getCentroid()) * facet.area();
}


std::unique_ptr<QuadratureRule> Tetrahedralization::clone() const
{
	return std::unique_ptr<QuadratureRule>(new Tetrahedralization(*this));
}

Real Tetrahedralization::apply(const Cell & cell, const std::function<Real(Point3D)> & integrand) const
{
	Real result = 0;
	auto & M_vertexIds = cell.getVerticesIds();
    assert(M_vertexIds.size() >= 4);

    if(M_vertexIds.size() > 4)
    {
		std::vector<Point3D> nodes;
        std::vector< std::vector<UInt> > facets;
        std::map<UInt, UInt> globalToLocal;

        // Count the number of vertices and facets that
        // define the approximation of the non-planar facets polyhedron
        const UInt nNodes = cell.verticesNumber();
        UInt addedNodes = 0;
        UInt nFacets = 0;
        auto & M_facetIds = cell.getFacetsIds();
        for(auto id : M_facetIds)
        {
            const UInt nodesFacet = cell.getMesh()->getFacetsVector()[id].verticesNumber();

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
        std::vector<UInt>::iterator it;

        nodes.resize( nNodes + addedNodes );
        facets.resize( nFacets );

        Point3D cellCenter(0., 0., 0.);

        for(i=0; i < nNodes; ++i)
        {
            nodes[i] = cell.getMesh()->getNodesVector()[M_vertexIds[i]];
            globalToLocal.insert(std::make_pair(M_vertexIds[i], i));
            cellCenter += nodes[i];
        }
        cellCenter /= nNodes;

        for(auto it=M_facetIds.begin(); it != M_facetIds.end(); ++it)
        {
			i = 0;
            auto& facet = cell.getMesh()->getFacetsVector()[*it];
            nodesFacet = facet.verticesNumber();

            if( nodesFacet > 3 )
            {
                Point3D center(0., 0., 0.);

                for(auto id : facet.getVerticesIds())
                {
                    center += cell.getMesh()->getNodesVector()[id];
                }
                center /= nodesFacet;
                nodes[nNodes + intNodesCount] = center;

                // loop over the triangles (except the last one) that settle the current facet
                for(j=0; j < nodesFacet - 1; ++j)
                {
                    // add the j-th triangle
                    facets[intFacetsCount].resize(3);
                    facets[intFacetsCount][0] = nNodes + intNodesCount;
                    facets[intFacetsCount][1] = globalToLocal[ facet.getVerticesIds()[j] ];
                    facets[intFacetsCount][2] = globalToLocal[ facet.getVerticesIds()[j+1] ];
                    intFacetsCount++;
                }
                // the last triangle that settles the current facet
                facets[intFacetsCount].resize(3);
                facets[intFacetsCount][0] = nNodes + intNodesCount;
                facets[intFacetsCount][1] = globalToLocal[ facet.getVerticesIds()[nodesFacet - 1] ];
                facets[intFacetsCount][2] = globalToLocal[ facet.getVerticesIds()[0] ];
                intFacetsCount++;

                intNodesCount++;
            }
            else
            {
                facets[intFacetsCount].resize(3);
                for(j=0; j < 3; ++j)
                    facets[intFacetsCount][j] = globalToLocal[ facet.getVerticesIds()[j] ];
                intFacetsCount++;
            }
            ++i;
        }

        // For each (triangular) facet we build a tetrahedron by adding a
        // fourth point at the center of the cell
        for (std::vector< std::vector<UInt> >::const_iterator facetIt = facets.begin();
            facetIt != facets.end(); ++facetIt)
        {
            std::vector<Point3D> tetrahedronNodes(4);
            tetrahedronNodes[0] = nodes[(*facetIt)[0]];
            tetrahedronNodes[1] = nodes[(*facetIt)[1]];
            tetrahedronNodes[2] = nodes[(*facetIt)[2]];
            tetrahedronNodes[3] = cellCenter;
            Real partialVolume = FVCode3D::tetrahedronVolume(tetrahedronNodes);
            Point3D centrioidTetra = ( tetrahedronNodes[0] + tetrahedronNodes[1] + tetrahedronNodes[2] + tetrahedronNodes[3] ) / 4.;
            result += partialVolume*integrand(centrioidTetra);
        }
    }
    else
    {
        std::vector<Point3D> nodes(4);
        nodes[0] = cell.getMesh()->getNodesVector()[M_vertexIds[0]];
        nodes[1] = cell.getMesh()->getNodesVector()[M_vertexIds[1]];
        nodes[2] = cell.getMesh()->getNodesVector()[M_vertexIds[2]];
        nodes[3] = cell.getMesh()->getNodesVector()[M_vertexIds[3]];

        Point3D centroid = ( nodes[0] + nodes[1] + nodes[2] + nodes[3] ) / 4.;
        result = cell.getVolume()*integrand(centroid);
    }
    return result;
}

Real Tetrahedralization::apply(const Facet & facet, const std::function<Real(Point3D)> & integrand) const
{
	Real result = 0.;
	auto & M_vertexIds = facet.getVerticesIds();
	assert(M_vertexIds.size() >= 3);
	
    if(M_vertexIds.size() > 3)
    {
        Point3D center(0., 0., 0.);
        for(auto id : M_vertexIds)
        {
            center += facet.getMesh()->getNodesVector()[id];
        }
        center /= M_vertexIds.size();

        auto it1 = M_vertexIds.begin();
        auto it2 = it1;
        it2++;
        Real area = 0.;

        for( ; it2 != M_vertexIds.end(); ++it1, ++it2)
        {
            const Point3D B = facet.getMesh()->getNodesVector()[*it1];
            const Point3D C = facet.getMesh()->getNodesVector()[*it2];

            area = FVCode3D::triangleArea(center, B, C);
            result += area * integrand( FVCode3D::triangleCentroid(center,B,C) );
        }

        const Point3D B = facet.getMesh()->getNodesVector()[ M_vertexIds[M_vertexIds.size()-1] ];
        const Point3D C = facet.getMesh()->getNodesVector()[ M_vertexIds[0] ];

        area = FVCode3D::triangleArea(center, B, C);
        result += area * integrand( FVCode3D::triangleCentroid(center,B,C) );
        
    }
    else
    {
        const Point3D A(facet.getMesh()->getNodesVector()[M_vertexIds[0]]);
        const Point3D B(facet.getMesh()->getNodesVector()[M_vertexIds[1]]);
        const Point3D C(facet.getMesh()->getNodesVector()[M_vertexIds[2]]);
        Real area = FVCode3D::triangleArea(A, B, C);
        result = area * integrand( FVCode3D::triangleCentroid(A,B,C) );
    }	
    return result;
}

} // namespace FVCode3D
