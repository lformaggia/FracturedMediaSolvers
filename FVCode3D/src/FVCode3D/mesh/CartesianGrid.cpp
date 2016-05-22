 /*!
 * @file CartesianGrid.cpp
 * @brief Class that generate a hexahedral structured (Cartesian) grid (definitions).
 */

#include <FVCode3D/mesh/Mesh3D.hpp>
#include <FVCode3D/property/Properties.hpp>
#include <FVCode3D/mesh/CartesianGrid.hpp>
#include <FVCode3D/mesh/MeshUtility.hpp>

namespace FVCode3D
{

void CartesianGrid::
generate( bool fracturesOn,
          const std::function<Real(const Point3D&)> & _surface )
{
    const Real Lx = M_data->getLx();
    const Real Ly = M_data->getLy();
//    const Real Lz = M_data->getLz();
    const UInt Nx = M_data->getNx();
    const UInt Ny = M_data->getNy();
    const UInt Nz = M_data->getNz();
    const Real Sx = M_data->getSx();
    const Real Sy = M_data->getSy();
    const Real Sz = M_data->getSz();
    const Real Rz = M_data->getRz() / 180. * _PI_;
    const bool noise = M_data->noiseOn();
    const Data::NoiseOn noiseOn = M_data->getNoiseOn();
    const Real mean = M_data->getMeanNormalDistribution();
    const Real stDev = M_data->getStDevNormalDistribution(); // 5.0e-2  7.5e-3 , 3.75e-3 , 1.875e-3 , 9.375e-4

    Real hx = Lx/Nx;
    Real hy = Ly/Ny;
//    Real hz = Lz/Nz;
    UInt nNodes = (Nx+1)*(Ny+1)*(Nz+1);
    UInt i,j,k;

    UInt bcId, zone, maxZone = 0, count = 0;

    std::vector<UInt> tmp, tmpFacets;

    std::vector<Point3D> & nodesRef = M_mesh.getNodesVector();
    std::map<UInt, Mesh3D::Facet3D> & facetsRef = M_mesh.getFacetsMap();
    std::map<UInt, Mesh3D::Cell3D> & cellsRef = M_mesh.getCellsMap();

    nodesRef.reserve(nNodes);

    // create nodes
    for(k=0; k <= Nz; ++k)
    {
        for(j=0; j <= Ny; ++j)
        {
            for(i=0; i <= Nx; ++i)
            {
                const Real abscissa = hx*i*cos(Rz) - hy*j*sin(Rz) + Sx;
                const Real ordinate = hx*i*sin(Rz) + hy*j*cos(Rz) + Sy;
                const Real Lz = _surface( Point3D( abscissa, ordinate, 0. ) );
                const Real quota = Lz/Nz*k + Sz;
                nodesRef.emplace_back( abscissa, ordinate, quota ); // Point3D
            }
        }
    }

    if(noise && (noiseOn == Data::NoiseOn::All))
    {
        addNoiseToPoint(M_mesh, mean, stDev);
    }

    tmp.resize(4);

    // create facets parallel to x axis
    for(i=0; i <= Nx; ++i)
    {
        for(k=0; k < Nz; ++k)
        {
            for(j=0; j < Ny; ++j)
            {
                tmp[0] = i + (Nx+1) * j + (Nx+1)*(Ny+1) * k;
                tmp[1] = i + (Nx+1) * j + (Nx+1)*(Ny+1) * (k+1);
                tmp[2] = i + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * (k+1);
                tmp[3] = i + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * k;

                bcId = ((i == 0) || (i == Nx)) ? 1 : 0;
                zone = 0;

                facetsRef.emplace( std::piecewise_construct, std::forward_as_tuple(count), std::forward_as_tuple(&M_mesh, tmp, (zone)*static_cast<UInt>(fracturesOn), bcId) );
                ++count;
            }
        }
    }

    // create facets parallel to y axis
    for(j=0; j <= Ny; ++j)
    {
        for(i=0; i < Nx; ++i)
        {
            for(k=0; k < Nz; ++k)
            {
                tmp[0] = i + (Nx+1) * j + (Nx+1)*(Ny+1) * k;
                tmp[1] = i + (Nx+1) * j + (Nx+1)*(Ny+1) * (k+1);
                tmp[2] = i+1 + (Nx+1) * j + (Nx+1)*(Ny+1) * (k+1);
                tmp[3] = i+1 + (Nx+1) * j + (Nx+1)*(Ny+1) * k;

                bcId = (( j==0 ) || (j == Ny)) ? 1 : 0;
                zone = 0;

                facetsRef.emplace( std::piecewise_construct, std::forward_as_tuple(count), std::forward_as_tuple(&M_mesh, tmp, (zone)*static_cast<UInt>(fracturesOn), bcId) );
                ++count;
            }
        }
    }

    // create facets parallel to z axis
    for(k=0; k <= Nz; ++k)
    {
        for(j=0; j < Ny; ++j)
        {
            for(i=0; i < Nx; ++i)
            {
                tmp[0] = i + (Nx+1) * j + (Nx+1)*(Ny+1) * k;
                tmp[1] = i+1 + (Nx+1) * j + (Nx+1)*(Ny+1) * k;
                tmp[2] = i+1 + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * k;
                tmp[3] = i + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * k;

                bcId = ((k == 0) || (k == Nz)) ? 1 : 0;
                zone = 0;

                facetsRef.emplace( std::piecewise_construct, std::forward_as_tuple(count), std::forward_as_tuple(&M_mesh, tmp, (zone)*static_cast<UInt>(fracturesOn), bcId) );
                ++count;
            }
        }
    }

    tmpFacets.resize(6);

    count = 0;

    M_mesh.buildNodesToFacetMap();

    // create cells
    for(k=0; k < Nz; ++k)
    {
        for(j=0; j < Ny; ++j)
        {
            for(i=0; i < Nx; ++i)
            {
                // bottom
                tmp[0] = i + (Nx+1) * j + (Nx+1)*(Ny+1) * k;
                tmp[1] = i+1 + (Nx+1) * j + (Nx+1)*(Ny+1) * k;
                tmp[2] = i+1 + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * k;
                tmp[3] = i + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * k;

                tmpFacets[0] = M_mesh.getFacetFromNodes(tmp);

                // top
                tmp[0] = i + (Nx+1) * j + (Nx+1)*(Ny+1) * (k+1);
                tmp[1] = i+1 + (Nx+1) * j + (Nx+1)*(Ny+1) * (k+1);
                tmp[2] = i+1 + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * (k+1);
                tmp[3] = i + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * (k+1);

                tmpFacets[1] = M_mesh.getFacetFromNodes(tmp);

                // left
                tmp[0] = i + (Nx+1) * j + (Nx+1)*(Ny+1) * k;
                tmp[1] = i + (Nx+1) * j + (Nx+1)*(Ny+1) * (k+1);
                tmp[2] = i+1 + (Nx+1) * j + (Nx+1)*(Ny+1) * (k+1);
                tmp[3] = i+1 + (Nx+1) * j + (Nx+1)*(Ny+1) * k;

                tmpFacets[2] = M_mesh.getFacetFromNodes(tmp);

                // right
                tmp[0] = i + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * k;
                tmp[1] = i + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * (k+1);
                tmp[2] = i+1 + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * (k+1);
                tmp[3] = i+1 + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * k;

                tmpFacets[3] = M_mesh.getFacetFromNodes(tmp);

                // back
                tmp[0] = i + (Nx+1) * j + (Nx+1)*(Ny+1) * k;
                tmp[1] = i + (Nx+1) * j + (Nx+1)*(Ny+1) * (k+1);
                tmp[2] = i + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * (k+1);
                tmp[3] = i + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * k;

                tmpFacets[4] = M_mesh.getFacetFromNodes(tmp);

                // front
                tmp[0] = i+1 + (Nx+1) * j + (Nx+1)*(Ny+1) * k;
                tmp[1] = i+1 + (Nx+1) * j + (Nx+1)*(Ny+1) * (k+1);
                tmp[2] = i+1 + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * (k+1);
                tmp[3] = i+1 + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * k;

                tmpFacets[5] = M_mesh.getFacetFromNodes(tmp);

                cellsRef.emplace( std::piecewise_construct, std::forward_as_tuple(count), std::forward_as_tuple(&M_mesh, tmpFacets, maxZone+1) );

                ++count;
            }
        }
    }
}

void CartesianGrid::extractBC(const Real theta)
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

            if(normal[compMax]>0.)
                it->second.setBorderID(2*compMax+1);
            else
                it->second.setBorderID(2*compMax+2);
        }
    }
}

void CartesianGrid::addFractures(const std::map<UInt,UInt> & facetIdToZone)
{
    std::map<UInt, Mesh3D::Facet3D> & facetsRef = M_mesh.getFacetsMap();

    FractureNetwork3D FN(M_mesh);
    std::vector<Fracture3D> fracturesVector;
    std::map<UInt, Fracture3D> fracturesMap;
    std::map<UInt, Fracture3D>::iterator itF;
    Point3D normal, center, centerFace;
    const UInt Nx = M_data->getNx();
    const UInt Ny = M_data->getNy();
    const UInt Nz = M_data->getNz();
    const bool noise = M_data->noiseOn();
    const Data::NoiseOn noiseOn = M_data->getNoiseOn();
    const Real mean = M_data->getMeanNormalDistribution();
    const Real stDev = M_data->getStDevNormalDistribution();

    for(std::map<UInt,UInt>::const_iterator it = facetIdToZone.begin(); it!= facetIdToZone.end(); ++it)
        if(facetsRef.find(it->first) != facetsRef.end())
            facetsRef[it->first].setZoneCode(it->second);

    for(std::map<UInt, Mesh3D::Facet3D>::iterator it = facetsRef.begin(); it != facetsRef.end(); ++it)
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

    if(noise && (noiseOn != Data::NoiseOn::All))
    {
        const bool initVal = (noiseOn == Data::NoiseOn::Matrix) ? true : false;
        const bool setVal = !initVal;
        std::vector<Point3D> & nodesRef = M_mesh.getNodesVector();
        std::map<UInt, Mesh3D::Cell3D> & cellsRef = M_mesh.getCellsMap();
        std::vector<bool> nodesWithNoise( nodesRef.size(), initVal );
        std::array<UInt,2> Cx={{0,Nx}};
        std::array<UInt,2> Cy={{0,Ny}};
        std::array<UInt,2> Cz={{0,Nz}};

        for(std::map<UInt, Fracture3D>::iterator it = fracturesMap.begin(); it != fracturesMap.end(); ++it)
        {
            for(auto facetId : it->second.getFractureFacetsId() )
            {
                for(auto nodeId : facetsRef[facetId].getVerticesVector())
                {
                    nodesWithNoise[nodeId] = setVal;
                }
            }
        }

        for(auto i : Cx)
        {
            for(auto j : Cy)
            {
                for(UInt k=0; k <= Nz; ++k)
                {
                    const Real index = i + (Nx+1) * j + (Nx+1) * (Ny+1) * k;
                    nodesWithNoise[index] = false;
                }
            }
            for(auto k : Cz)
            {
                for(UInt j=0; j <= Ny; ++j)
                {
                    const Real index = i + (Nx+1) * j + (Nx+1) * (Ny+1) * k;
                    nodesWithNoise[index] = false;
                }
            }
        }
        for(auto j : Cy)
        {
            for(auto k : Cz)
            {
                for(UInt i=0; i <= Nx; ++i)
                {
                    const Real index = i + (Nx+1) * j + (Nx+1) * (Ny+1) * k;
                    nodesWithNoise[index] = false;
                }
            }
        }

        std::vector<bool> nodesWithNoiseX( nodesRef.size(), false );
        std::vector<bool> nodesWithNoiseY( nodesRef.size(), false );
        std::vector<bool> nodesWithNoiseZ( nodesRef.size(), false );

        for(auto i : Cx)
        {
            for(UInt j=1; j < Ny; ++j)
            {
                for(UInt k=1; k < Nz; ++k)
                {
                    const Real index = i + (Nx+1) * j + (Nx+1) * (Ny+1) * k;
                    nodesWithNoiseX[index] = nodesWithNoise[index];
                    nodesWithNoise[index] = false;
                }
            }
        }
        for(auto j : Cy)
        {
            for(UInt i=1; i < Nx; ++i)
            {
                for(UInt k=1; k < Nz; ++k)
                {
                    const Real index = i + (Nx+1) * j + (Nx+1) * (Ny+1) * k;
                    nodesWithNoiseY[index] = nodesWithNoise[index];
                    nodesWithNoise[index] = false;
                }
            }
        }
        for(auto k : Cz)
        {
            for(UInt i=1; i < Nx; ++i)
            {
                for(UInt j=1; j < Ny; ++j)
                {
                    const Real index = i + (Nx+1) * j + (Nx+1) * (Ny+1) * k;
                    nodesWithNoiseZ[index] = nodesWithNoise[index];
                    nodesWithNoise[index] = false;
                }
            }
        }

        addNoiseToPoint(M_mesh, nodesWithNoise, mean, stDev);
        addNoiseToPoint(M_mesh, nodesWithNoiseX, std::vector<UInt>{1,2}, mean, stDev);
        addNoiseToPoint(M_mesh, nodesWithNoiseY, std::vector<UInt>{0,2}, mean, stDev);
        addNoiseToPoint(M_mesh, nodesWithNoiseZ, std::vector<UInt>{0,1}, mean, stDev);

        for(std::map<UInt, Mesh3D::Facet3D>::iterator it = facetsRef.begin(); it != facetsRef.end(); ++it)
        {
            it->second.computeCentroidAndNormalAndArea();
        }

        for(std::map<UInt, Mesh3D::Cell3D>::iterator it = cellsRef.begin(); it != cellsRef.end(); ++it)
        {
            it->second.computeVolumeAndCentroid();
        }
    }
}

void CartesianGrid::addBCAndFractures(const std::map<UInt,UInt> & facetIdToZone, const Real theta)
{
    extractBC(theta);
    addFractures(facetIdToZone);
}

} // namespace FVCode3D
