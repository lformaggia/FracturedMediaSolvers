/*!
 *	@file converter.hpp
 *	@brief Methods to convert format files.
 */

#ifndef CONVERTER_HPP_
#define CONVERTER_HPP_

#include "core/TypeDefinition.hpp"

class Mesh3D;
class PropertiesMap;

//! Add the lacking data to the standard TPFA format
/*!
 * Add the BCs and the fracture network
 * @param mesh reference to a Geometry::Mesh3D
 * @param properties reference to a Geometry::PropertiesMap
 * @param theta rotation angle along z-axis. It is used only to compute the BC ids. Default = 0
 * @pre import grid file by means of the TPFA Standard format
 */
void readTPFAStandardAsTPFAWithBC(Geometry::Mesh3D & mesh, Geometry::PropertiesMap & propertie, const Real theta = 0.);

//! Generate the BC ids from a standard TPFA file format
/*!
 * Add the BCs ids to the boundary facets.
 * The id is set by considering the maximum component and the sign of the normal of a boundary facet.
 * @param mesh reference to a Geometry::Mesh3D
 * @param properties reference to a Geometry::PropertiesMap
 * @param theta rotation angle along z-axis. It is used only to compute the BC ids. Default = 0
 * @pre import grid file by means of the TPFA Standard format
 */
void extractBC(Geometry::Mesh3D & mesh, Geometry::PropertiesMap & properties, const Real theta = 0.);

//! Save a grid file as an improved TPFA format, by adding the lacking data
/*!
 * Add the BCs and the fracture network as save in grid file format.
 * @param filename name of the file
 * @param mesh reference to a Geometry::Mesh3D
 * @param properties reference to a Geometry::PropertiesMap
 * @pre import grid file by means of the TPFA Standard format
 */
void convertFromTPFAStandardToTPFAWithBC(const std::string filename, Geometry::Mesh3D & mesh, Geometry::PropertiesMap & properties);

#endif /* CONVERTER_HPP_ */
