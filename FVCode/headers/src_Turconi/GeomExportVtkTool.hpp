 /*!
 *	@file GeomExportVtkTool.hpp
 *	@brief Some function to Export in Vtk the elementary geometrical objects.
 *
 *	@author Luca Turconi <luca.turconi@moxoff.com>
 *
 */ 

#ifndef GEOMEXPORTVTKTOOL_HPP_
#define GEOMEXPORTVTKTOOL_HPP_

#include "TypeDefinition.hpp"

#include<vector>

namespace Geometry{
	
	bool exportSetOfPointsVtk( const std::vector<Geometry::Point2D> & set,
							   const std::string & fileName );
							   
	bool exportSegment2DVtk( const Geometry::Segment2D & seg,
							 const std::string & fileName );
	
} // namespace Geometry

#endif /* GEOMEXPORTVTKTOOL_HPP_ */
