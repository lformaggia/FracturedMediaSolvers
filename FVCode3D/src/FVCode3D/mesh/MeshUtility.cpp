/*!
 * @file MeshUtility.hpp
 * @brief This unit contains some utilities concerning meshes (definitions).
 */

#include <FVCode3D/mesh/MeshUtility.hpp>
#include <FVCode3D/mesh/Mesh3D.hpp>

namespace FVCode3D
{

void addNoiseToPoint(Mesh3D & mesh, const Real mean, const Real stDev)
{
    std::default_random_engine generator;
    std::normal_distribution<Real> distribution(mean, stDev);

    for(auto& node : mesh.getNodesVector())
    {
        node.x() += distribution(generator);
        node.z() += distribution(generator);
        node.y() += distribution(generator);
    }
} // addNoiseToPoint

void addNoiseToPoint(Mesh3D & mesh, const std::vector<bool> & nodesWithNoise, const Real mean,
    const Real stDev)
{
    std::default_random_engine generator;
    std::normal_distribution<Real> distribution(mean, stDev);
    std::vector<Point3D> & nodesRef = mesh.getNodesVector();
    const UInt NoN ( mesh.getNumberOfNodes() );

    for(auto i = 0; i<NoN; ++i)
    {
        if(nodesWithNoise[i])
        {
            nodesRef[i].x() += distribution(generator);
            nodesRef[i].y() += distribution(generator);
            nodesRef[i].z() += distribution(generator);
        }
    }
} // addNoiseToPoint

void addNoiseToPoint(Mesh3D & mesh, const std::vector<bool> & nodesWithNoise, const std::vector<UInt> & coords,
    const Real mean, const Real stDev)
{
    std::default_random_engine generator;
    std::normal_distribution<Real> distribution(mean, stDev);
    std::vector<Point3D> & nodesRef = mesh.getNodesVector();
    const UInt NoN ( mesh.getNumberOfNodes() );

    for(auto i = 0; i<NoN; ++i)
    {
        if(nodesWithNoise[i])
        {
            for(auto coord : coords)
            {
                nodesRef[i][coord] += distribution(generator);
            }
        }
    }
} // addNoiseToPoint

} // namespace FVCode3D
