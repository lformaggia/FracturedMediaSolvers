/*!
 * @file converter.hpp
 * @brief Methods to convert format files.
 */

#ifndef CONVERTER_HPP_
#define CONVERTER_HPP_

#include <FVCode3D/core/TypeDefinition.hpp>

namespace FVCode3D
{

class Mesh3D;
class PropertiesMap;

//! Save a file in solver format (.fvg)
/*!
 * Save mesh, properties and fractures.
 * @param filename name of the file
 * @param mesh reference to a Mesh3D
 * @param properties reference to a PropertiesMap
 * @pre import the file
 */
void saveAsSolverFormat(const std::string filename, Mesh3D & mesh, PropertiesMap & properties);

//! Save a file in Medit format (.mesh)
/*!
 * Save mesh and fractures(as a suitable label).
 * @param filename name of the file
 * @param mesh reference to a Mesh3D
 * @pre import the file
 */
void saveAsMeditFormat(const std::string filename, Mesh3D & mesh);

//! Save a file in OpenFOAM format (polyMesh)
/*!
 * Save mesh and fractures(as a suitable label).
 * @param filename name of the file
 * @param mesh reference to a Mesh3D
 * @pre import the file
 */
void saveAsOpenFOAMFormat(const std::string filename, Mesh3D & mesh);

} // namespace FVCode3D
#endif /* CONVERTER_HPP_ */
