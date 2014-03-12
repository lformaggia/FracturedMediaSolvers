/*!
 *	@file converter.hpp
 *	@brief Methods to convert format files.
 */

#ifndef CONVERTER_HPP_
#define CONVERTER_HPP_

#include "core/TypeDefinition.hpp"

class Mesh3D;
class PropertiesMap;

//! Save a file as a file for solver format
/*!
 * Save mesh, fractures and properties.
 * Further, add the BCs and the fracture network.
 * @param filename name of the file
 * @param mesh reference to a Geometry::Mesh3D
 * @param properties reference to a Geometry::PropertiesMap
 * @pre import the file
 */
void saveAsSolverFormat(const std::string filename, Geometry::Mesh3D & mesh, Geometry::PropertiesMap & properties);

#endif /* CONVERTER_HPP_ */
