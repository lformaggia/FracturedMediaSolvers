 /*!
 *	@file Dimension.hpp
 *	@brief Class for type trait.
 *
 *	@author Francesco Della Porta 
 *
 */ 
#ifndef DIMENSION_HPP_
#define DIMENSION_HPP_
 
#include "TypeDefinition3D.hpp"

#include<set>
#include<vector>
#include<utility>
#include<map>



namespace Geometry{
	
/*!
	@class Dimension
	This class is used as a template parameter in order to decide at compilation time if a mesh is bi- or tri-dimensional. This class is a template class and admits as a template parameter just 2 or 3 for, respectively, 2D or 3D meshes.
*/
template<int N>
struct Dimension
{};

/*!
	@class Dimension<2>
	This class is the specialization of the class Dimension for the parameter 2. It has in it many typedefinitions usefull to distinguish 2D from 3D objects*/
template<>
struct Dimension<2>
{
	/*!
		@typedef Generic_Point
   		This typedefinition permits to handle Point2D as a Generic_Point.
   	*/
	typedef Geometry::Point2D Generic_Point;
	/*!
		@typedef Generic_Vector
   		This typedefinition permits to handle a Vector2D as a Generic_Vector.
   	*/
	typedef Geometry::Vector2D Generic_Vector;
	/*!
		@typedef Generic_Segment
   		This typedefinition permits to handle a Segment2D as a Generic_Segment.
   	*/
	typedef Geometry::Segment2D Generic_Segment;

	/*!
		@typedef Generic_Facet
   		This typedefinition permits to handle a Mesh2D::Edge2D as a Generic_Facet.
   	*/
	typedef Geometry::Mesh2D::Edge2D Generic_Facet;
	/*!
		@typedef Generic_Cell
   		This typedefinition permits to handle a Mesh2D::Cell2D as a Generic_Cell.
   	*/
	typedef Geometry::Mesh2D::Cell2D Generic_Cell;
	/*!
		@typedef Generic_Mesh
   		This typedefinition permits to handle a Mesh2D as a Generic_Mesh.
   	*/
	typedef Geometry::Mesh2D Generic_Mesh;

	/*!
		@typedef FractureNetwork
   		This typedefinition permits to handle a FractureNetwork2D as a generic FractureNetwork.
   	*/
	typedef Geometry::FractureNetwork2D FractureNetwork;
	/*!
		@typedef Fracture_Juncture
   		This typedefinition permits to handle a Point-Id as a juncture between to facets.
   	*/
	typedef UInt Fracture_Juncture;
	/*!
		@typedef Generic_Border
   		This typedefinition permits to handle a Segment2D as a Generic_Border of the domain.
   	*/
	typedef Geometry::Segment2D Generic_Border;
	enum {dim=2};
};

/*!
	@class Dimension<3>
	This class is the specialization of the class Dimension for the parameter 3. It has in it many typedefinitions usefull to distinguish 2D from 3D objects*/
template<>
struct Dimension<3>
{
	/*!
		@typedef Generic_Point
   		This typedefinition permits to handle Point3D as a Generic_Point.
   	*/
	typedef Geometry::Point3D Generic_Point;
	/*!
		@typedef Generic_Vector
   		This typedefinition permits to handle a Vector3D as a Generic_Vector.
   	*/
	typedef Geometry::Vector3D Generic_Vector;
	/*!
		@typedef Generic_Segment
   		This typedefinition permits to handle a Segment3D as a Generic_Segment.
   	*/
	typedef Geometry::Segment3D Generic_Segment;

	/*!
		@typedef Generic_Facet
   		This typedefinition permits to handle a Mesh3D::Facet3D as a Generic_Facet.
   	*/
	typedef Geometry::Mesh3D::Facet3D Generic_Facet;
	/*!
		@typedef Generic_Cell
   		This typedefinition permits to handle a Mesh3D::Cell3D as a Generic_Cell.
   	*/
	typedef Geometry::Mesh3D::Cell3D Generic_Cell;
	/*!
		@typedef Generic_Mesh
   		This typedefinition permits to handle a Mesh3D as a Generic_Mesh.
   	*/
	typedef Geometry::Mesh3D Generic_Mesh;

	/*!
		@typedef FractureNetwork
   		This typedefinition permits to handle a FractureNetwork3D as a generic FractureNetwork.
   	*/
	typedef std::pair<UInt,UInt> Fracture_Juncture;
	/*!
		@typedef Fracture_Juncture
   		This typedefinition permits to handle a pair of Points as a juncture between to facets.
   	*/
	typedef Geometry::FractureNetwork3D FractureNetwork;
	/*!
		@typedef Generic_Border
   		This typedefinition permits to handle a vector of Point3D as a Generic_Border of the domain.
   	*/
	typedef std::vector<Geometry::Point3D> Generic_Border;
	enum {dim=3};
};

}
#endif
