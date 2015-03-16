/*!
 * @file MeshUtility.hpp
 * @brief This unit contains some utilities concerning meshes.
 */

#ifndef MESHUTILITY_HPP_
#define MESHUTILITY_HPP_

#include <FVCode3D/core/TypeDefinition.hpp>
#include <FVCode3D/core/Data.hpp>

namespace FVCode3D
{

class Mesh3D;

//! Add noise to all points
/*!
 * Add noise to the points following a normal distribution with mean @a mean and standard deviation @a stDev
 * @param mesh reference to a Mesh3D
 * @param mean mean. Default = 0.
 * @param stDev standard deviation. Default = 1.
 */
void addNoiseToPoint(Mesh3D & mesh, const Real mean = 0., const Real stDev = 1.);

//! Add noise to a set of points
/*!
 * Add noise to the points following a normal distribution with mean @a mean and standard deviation @a stDev
 * @param mesh reference to a Mesh3D
 * @param nodesWithNoise vector that indicates on which nodes the noise is applied
 * @param mean mean. Default = 0.
 * @param stDev standard deviation. Default = 1.
 */
void addNoiseToPoint(Mesh3D & mesh, const std::vector<bool> & nodesWithNoise, const Real mean = 0.,
    const Real stDev = 1.);

//! Add noise to a set of points
/*!
 * Add noise to the points following a normal distribution with mean @a mean and standard deviation @a stDev
 * @param mesh reference to a Mesh3D
 * @param nodesWithNoise vector that indicates on which nodes the noise is applied
 * @param coords coordinates on which the noise is applied
 * @param mean mean. Default = 0.
 * @param stDev standard deviation. Default = 1.
 */
void addNoiseToPoint(Mesh3D & mesh, const std::vector<bool> & nodesWithNoise, const std::vector<UInt> & coords,
    const Real mean = 0., const Real stDev = 1.);

} // namespace FVCode3D

#endif // MESHUTILITY_HPP_
