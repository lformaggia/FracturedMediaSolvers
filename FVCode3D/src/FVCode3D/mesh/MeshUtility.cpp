/*!
 * @file MeshUtility.hpp
 * @brief This unit contains some utilities concerning meshes (definitions).
 */

#include <FVCode3D/mesh/MeshUtility.hpp>
#include <FVCode3D/mesh/Mesh3D.hpp>

// Set this macro to apply the noise to all nodes...
#define NOISE_TO_ALL_NODES
// ... otherwise it applies only to internal nodes
// WARNING: only if the mesh is a Cartesian Grid.
#ifndef NOISE_TO_ALL_NODES
#define NOISE_TO_INTERNAL_NODES
#endif // NOISE_TO_ALL_NODES

namespace FVCode3D
{

void addNoiseToPoint(Mesh3D & mesh, const Real mean, const Real stDev)
{
    std::default_random_engine generator;
    std::normal_distribution<Real> distribution(mean, stDev);

#ifdef NOISE_TO_ALL_NODES
    for(auto& node : mesh.getNodesVector())
    {
        node.x() += distribution(generator);
        node.z() += distribution(generator);
        node.y() += distribution(generator);
    }
#endif // NOISE_TO_ALL_NODES

#ifdef NOISE_TO_INTERNAL_NODES
    std::vector<Point3D> & nodesRef = mesh.getNodesVector();
    for(j=1; j < Ny; ++j)
    {
        for(i=1; i < Nx; ++i)
        {
            for(k=1; k < Nz; ++k)
            {
                const Real varX = distribution(generator);
                const Real varY = distribution(generator);
                const Real varZ = distribution(generator);
                const Real index = i + (Nx+1) * j + (Nx+1) * (Ny+1) * k;
                nodesRef[index].x() += varX;
                nodesRef[index].y() += varY;
                nodesRef[index].z() += varZ;
            }
        }
    }

    std::array<UInt,2> Cx={{0,Nx}};
    std::array<UInt,2> Cy={{0,Ny}};
    std::array<UInt,2> Cz={{0,Nz}};

    for(auto j : Cy)
    {
        for(i=1; i < Nx; ++i)
        {
            for(k=1; k < Nz; ++k)
            {
                const Real varX = distribution(generator);
                const Real varZ = distribution(generator);
                const Real index = i + (Nx+1) * j + (Nx+1) * (Ny+1) * k;
                nodesRef[index].x() += varX;
                nodesRef[index].z() += varZ;
            }
        }
    }

    for(auto i : Cx)
    {
        for(j=1; j < Ny; ++j)
        {
            for(k=1; k < Nz; ++k)
            {
                const Real varY = distribution(generator);
                const Real varZ = distribution(generator);
                const Real index = i + (Nx+1) * j + (Nx+1) * (Ny+1) * k;
                nodesRef[index].y() += varY;
                nodesRef[index].z() += varZ;
            }
        }
    }

    for(auto k : Cz)
    {
        for(i=1; i < Nx; ++i)
        {
            for(j=1; j < Ny; ++j)
            {
                const Real varX = distribution(generator);
                const Real varY = distribution(generator);
                const Real index = i + (Nx+1) * j + (Nx+1) * (Ny+1) * k;
                nodesRef[index].x() += varX;
                nodesRef[index].y() += varY;
            }
        }
    }
#endif // NOISE_TO_INTERNAL_NODES
} // addNoiseToPoint

} // namespace FVCode3D
