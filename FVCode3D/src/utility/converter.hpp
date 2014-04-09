/*!
 *	@file converter.hpp
 *	@brief Methods to convert format files.
 */

#ifndef CONVERTER_HPP_
#define CONVERTER_HPP_

#include "core/TypeDefinition.hpp"

class Mesh3D;
class PropertiesMap;

//! Save a file in solver format (.fvg)
/*!
 * Save mesh, properties and fractures.
 * @param filename name of the file
 * @param mesh reference to a Geometry::Mesh3D
 * @param properties reference to a Geometry::PropertiesMap
 * @pre import the file
 */
void saveAsSolverFormat(const std::string filename, Geometry::Mesh3D & mesh, Geometry::PropertiesMap & properties);

//! Save a file in Medit format (.mesh)
/*!
 * Save mesh and fractures(as a suitable label).
 * @param filename name of the file
 * @param mesh reference to a Geometry::Mesh3D
 * @pre import the file
 */
void saveAsMeditFormat(const std::string filename, Geometry::Mesh3D & mesh);

#endif /* CONVERTER_HPP_ */
