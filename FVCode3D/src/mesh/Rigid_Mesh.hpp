 /*!
 *	@file Rigid_Mesh.hpp
 *	@brief Class for unstructured mesh.
 */

#ifndef RIGID_MESH_HPP_
#define RIGID_MESH_HPP_
 
#include "core/TypeDefinition.hpp"
#include "mesh/Mesh3D.hpp"

#include <cmath>
#include <set>
#include <vector>
#include <tuple>
#include <algorithm>
#include <map>
#include <utility>

namespace Geometry{

class PropertiesMap;

//! Class that allows to handle a Mesh3D for a differential problem
/*!
	@class Rigid_Mesh
   	This class implements a container for a Geometry::Mesh3D.
   	The Rigid_Mesh container adapts the Geometry::Mesh3D in order to be efficient for the discretization of differential problems.
*/
class Rigid_Mesh{
public:

	//! Typedef for Geometry::Point3D
	/*!
		@typedef Generic_Point
   		This type definition permits to treat Geometry::Point3D as a Generic_Point.
   	*/
	typedef Point3D Generic_Point;
	//! Typedef for Geometry::Point3D
	/*!
		@typedef Generic_Vector
   		This type definition permits to treat Geometry::Point3D as a Generic_Vector.
   	*/
	typedef Point3D Generic_Vector;
	//! Typedef for Geometry::Mesh3D
	/*!
		@typedef Generic_Mesh
   		This type definition permits to treat Mesh3D as a Generic_Mesh.
   	*/
	typedef Mesh3D Generic_Mesh;
	//! Typedef for Geometry::Mesh3D::Facet3D
	/*!
		@typedef Generic_Facet
   		This type definition permits to treat Mesh3D::Facet3D as a Generic_Facet.
   	*/
	typedef Mesh3D::Facet3D Generic_Facet;
	//! Typedef for Geometry::Mesh3D::Cell3D
	/*!
		@typedef Generic_Cell
   		This type definition permits to treat Mesh3D::Cell3D as a Generic_Cell.
   	*/
	typedef Mesh3D::Cell3D Generic_Cell;
	//! Typedef for std::pair<UInt,UInt>
	/*!
		@typedef Generic_Edge
   		This type definition permits to treat a pair of Points as a Generic_Edge.
   	*/
	typedef std::pair<UInt,UInt> Generic_Edge;
	//! Typedef for std::pair<UInt,UInt>
	/*!
		@typedef Fracture_Juncture
   		This type definition permits to handle a pair of Points as a juncture between to facets.
   	*/
	typedef std::pair<UInt,UInt> Fracture_Juncture;
	//! Typedef for std::pair<UInt,UInt>
	/*!
		@typedef Fracture_Tip
   		This type definition permits to handle a pair of Points as a tip of a facet.
   	*/
	typedef std::pair<UInt,UInt> Fracture_Tip;
	//! Typedef for std::vector<Point3D>
	/*!
		@typedef Generic_Border
   		This type definition permits to handle a vector of Point3D as a Generic_Border of the domain.
   	*/
	typedef std::vector<Point3D> Generic_Border;
	//! Typedef for std::map<UInt,Mesh3D::Facet3D>
	/*!
		@typedef Generic_FacetsContainer
	   	This type definition permits to treat the container of the Facets of a Mesh3D as Generic_FacetsContainer.
	 */
	typedef std::map<UInt,Mesh3D::Facet3D> Generic_FacetsContainer;
	//! Typedef for Generic_Point
	/*!
		@typedef Point
   		This type definition permits to treat Generic_Point as a Point.
   	*/
	typedef Generic_Point Point;

	//! Class that implements a Edge
	/*!
		@class Edge
   		This class implements the concept of Edge of a Rigid_Mesh.
   	*/
	class Edge{
	public:

		//! @name Constructor & Destructor
		//@{

		//! Constructor for an Edge from a pair of points.
		/*!
			@param generic_edge is a edge
			@param mesh is a constant pointer to the mesh to which the Edge belongs
			@param id is the Id of the edge in the Rigid_Mesh
		*/
		Edge (const Generic_Edge & generic_edge, Geometry::Rigid_Mesh * const mesh, const UInt id);

		//! Copy constructor for a Edge from a Edge of another Rigid_Mesh
		/*!
		 * @param edge reference to a Edge
		 * @param mesh is a constant pointer to the mesh to which the Edge belongs
		 */
		Edge (const Edge & edge, Geometry::Rigid_Mesh * const mesh);

		//! Copy constructor for a Edge from a Edge
		/*!
		 * @param edge reference to a Edge
		 */
		Edge (const Edge & edge);

		//! No empty constructor
		Edge () = delete;

		//! Destructor
		~Edge () = default;

		//@}

		//! @name Get Methods
		//@{

		//! Get vertexes Ids (const)
		/*!
		 * @return a reference to a pair of the edge vertexes ids
		 */
		const Generic_Edge & getEdge () const
			{ return M_edge; }

		//! Get mesh
		/*!
		 * @return a const pointer to the related mesh
		 */
		Geometry::Rigid_Mesh * getMesh () const
			{ return M_mesh; }

		//! Get Id (const)
		/*!
		 * @return the id of the Edge in the containing Rigid_Mesh
		 */
		UInt getId() const
			{ return M_id; }
		//@}

		//! Get separated Facets ids (const)
		/*!
		 * @return a vector reference with the ids of the facets separated from the Edge
		 */
		const std::vector<UInt> & getSeparatedFacetsIds () const
			{ return M_separatedFacetsIds; }

		//! Get centroid (const)
		/*!
		 * @return reference to the centroid of the Edge
		 */
		const Generic_Point getCentroid () const;

		//! Get the length of the edge (const)
		/*!
		 * @return the length of the Edge
		 */
		Real length() const;

		//! @name Methods
		//@{
		//! Display general information about the Edge
		/*!
		 * @param out specify the output format (std::cout by default)
		 */
		void showMe (std::ostream & out=std::cout) const;

		//@}

		friend Rigid_Mesh;

	protected:
		//! The pair of the ids of the Edge's Vertexes
		const Generic_Edge M_edge;
		//! The pointer to the Rigid_Mesh containing the Edge
		Geometry::Rigid_Mesh * M_mesh;
		//! The Edge id in the containing Rigid_Mesh
		UInt M_id;
		//! The vector of the ids of the Facets separated by Edge
		std::vector<UInt> M_separatedFacetsIds;

	};

	//! Class that implements a Facet
	/*!
		@class Facet
   		This class implements the concept of Facet of a Rigid_Mesh.
   	*/
	class Facet{

	public:

		//! @name Constructor & Destructor
		//@{

		//! Constructor for a Facet from a facet of a Mesh3D.
		/*!
			@param generic_facet is a facet of a Mesh3D.
			@param mesh is a constant pointer to the mesh to which the Facet belongs
			@param old_to_new_map is a map from the old Id of the Neighbour-Cells in the Mesh3D to the new Id of those cells in the Rigid_Mesh
			@param m_id is the Id of the facet in the Rigid_Mesh
		*/
		Facet (const Generic_Facet & generic_facet, Geometry::Rigid_Mesh * const mesh, const std::map<UInt, UInt> &old_to_new_map, const UInt m_id);

		//! Copy constructor for a Facet from a Facet of another Rigid_Mesh
		/*!
		 * @param facet reference to a Facet
		 * @param mesh is a constant pointer to the mesh to which the Facet belongs
		 */
		Facet (const Facet& facet, Geometry::Rigid_Mesh * const mesh);

		//! Copy constructor for a Facet from a Facet
		/*!
		 * @param facet reference to a Facet
		 */
		Facet (const Facet& facet);

		//! No empty constructor
		Facet () = delete;

		//! Destructor
		~Facet () = default;

		//@}

		//! @name Get Methods
		//@{
		
		//! Get mesh
		/*!
		 * @return a const pointer to the related mesh
		 */
		Geometry::Rigid_Mesh * getMesh () const
			{ return M_mesh; }

		//! Get vertexes Ids (const)
		/*!
		 * @return a reference to a vector of the facet vertexes ids
		 */
		const std::vector<UInt> & getVertexesIds () const
			{ return M_vertexesIds; }

		//! Get Vertexes number (const)
		/*!
		 * @return the number of the vertexes of the Facet
		 */
		UInt vertexesNumber() const
			{ return M_vertexesIds.size (); }

		//! Get centroid (const)
		/*!
		 * @return reference to the centroid of the Facet
		 */
		const Generic_Point & getCentroid () const
			{ return M_centroid; }

		//! Get unsigned Normal (const)
		/*!
		 * @return the unsigned normal vector to the Facet
		 */
		Generic_Vector getUnsignedNormal () const
			{ return M_unsignedNormal; }

		//! Get separated Cells ids (const)
		/*!
		 * @return a vector reference with the ids of the cells separated from the Facet
		 */
		const std::vector<UInt> & getSeparatedCellsIds () const
			{ return M_separatedCellsIds; }

		//! Get the area of the facet (const)
		/*!
		 * @return the area of the Facet
		 */
		Real area() const
			{ return M_area; }

		//! Get Id (const)
		/*!
		 * @return the id of the Facet in the containing Rigid_Mesh
		 */
		UInt getId() const
			{ return M_id; }
		//@}

		//! @name Methods
		//@{
		//! Display general information about the Facet
		/*!
		 * @param out specify the output format (std::cout by default)
		 */	
		void showMe (std::ostream & out=std::cout) const;

		//@}

	protected:
		//! The pointer to the Rigid_Mesh containing the Facet
		Geometry::Rigid_Mesh * M_mesh;
		//! The Facet id in the containing Rigid_Mesh
		UInt M_id;
		//! The vector of the ids of the Facet's Vertexes
		std::vector<UInt> M_vertexesIds;
		//! The vector of the ids of the Cells separated by Facet
		std::vector<UInt> M_separatedCellsIds;
		//! The area of the Facet
		Real M_area;
		//! The centroid of the Facet
		Generic_Point M_centroid;
		//! The unsigned normal-vector to the Facet
		Generic_Vector M_unsignedNormal;

	};

	//! Class that implements a Cell
	/*!
		@class Cell
   		This class implements the concept of Cell of a Rigid_Mesh.
   	*/
	class Cell{

	public:
		//! @name Constructor & Destructor
		//@{

		//! Constructor for a Cell from a cell of a Mesh2D or a Mesh3D.
		/*!
			@param generic_cell is a facet of a Mesh3D.
			@param mesh is a constant pointer to the mesh to which the Cell belongs
			@param m_id is the id of the cell in the Rigid_Mesh
		*/
		Cell (const Generic_Cell & generic_cell, Geometry::Rigid_Mesh * const mesh, const UInt m_id);

		//! Copy constructor for a Cell from a Cell of another Rigid_Mesh
		/*!
		 * @param cell reference to a Cell
		 * @param mesh is a constant pointer to the mesh to which the Cell belongs
		 */
		Cell (const Cell & cell, Geometry::Rigid_Mesh * const mesh);

		//! Copy constructor for a Cell from a Cell
		/*!
		 * @param cell reference to a Cell
		 */
		Cell (const Cell & cell);

		//! No empty constructor
		Cell () = delete;

		//! Destructor
		~Cell() = default;

		//@}

		//! @name Get Methods
		//@{
		
		//! Get mesh
		/*!
		 * @return a const pointer to the mesh
		 */
		Geometry::Rigid_Mesh * getMesh () const
			{ return M_mesh; }

		//! Zone Code (const)
		/*!
		 * @return the zone code of the represented fractures
		 */
		UInt getZoneCode () const
			{return M_zoneCode;}

		//! Get Vertexes ids (const)
		/*!
		 * @return a reference to a vector which contains the ids of the Cell's vertices
		 */
		const std::vector<UInt> & getVertexesIds () const
			{ return M_vertexesIds; }

		//! Get Facets ids (const)
		/*!
		 * @return a reference to a vector which contains the ids of the Facets
		 */
		const std::vector<UInt> & getFacetsIds () const
			{ return M_facetsIds; }

		//! Get Neighbors ids (const)
		/*!
		 * @return a reference to a vector which contains the ids of the neighour Cells
		 */
		const std::vector<UInt> & getNeighborsIds () const
			{ return M_neighborsIds; }

		//! Get centroid (const)
		/*!
		 * @return a reference to the centroid of the Cell
		 */
		const Generic_Point & getCentroid () const
			{ return M_centroid; }

		//! Get Vertexes number (const)
		/*!
		 * @return the number of the vertexes of the Cell
		 */
		UInt vertexesNumber() const
			{ return M_vertexesIds.size (); }

		//! Get Facets number (const)
		/*!
		 * @return the number of the facets of the Cell
		 */
		UInt facetsNumber() const
			{ return M_facetsIds.size (); }

		//! Get neighbours number (const)
		/*!
		 * @return the number of neighbours of the Cell
		 */
		UInt neighborsNumber() const
			{ return M_neighborsIds.size (); }

		//! Get volume (const)
		/*!
		 * @return the volume of the Cell
		 */
		Real getVolume() const
			{ return M_volume; }

		//! Get cell id (const)
		/*!
		 * @return The id of the Cell in the containing Rigid_Mesh
		 */
		UInt getId() const
			{ return M_Id; }
		//@}

		//! @name Methods
		//@{
		//! Display general information about the Cell
		/*!
		 * @param out specify the output format (std::cout by default)
		 */	
		void showMe (std::ostream & out=std::cout) const;

		//@}

		friend class Rigid_Mesh;

	protected:
		//! The pointer to the Rigid_Mesh containing the Cell
		Geometry::Rigid_Mesh * M_mesh;
		//! The Cell Id in the containing Rigid_Mesh
		UInt M_Id;
		//! Zone code of the Cell
		UInt M_zoneCode;
		//! The vector of the ids of the Cell's Vertexes
		std::vector<UInt> M_vertexesIds;
		//! The vector of the ids of the Facets
		std::vector<UInt> M_facetsIds;
		//! The vector of the ids of the Cells neighbors
		std::vector<UInt> M_neighborsIds;
		//! The centroid of the Cell
		Generic_Point M_centroid;
		//! The volume of the Cell
		Real M_volume;
	};

	//! Class that handles the several types of edges
	/*!
		@class Edge_ID
   		This class is a container of a Edge-Id.
		This class is a base class.
   	*/
	class Edge_ID{

	protected:
		//! The Id of the referred Edge
		UInt Edge_Id;
		//! The pointer to the Rigid_Mesh containing the Edge
		Rigid_Mesh * M_mesh;

	public:
		//! @name Constructor & Destructor
		//@{

		//! Constructor for a Edge_ID.
		/*!
			@param edge_id is the id of a Edge of a Rigid_Mesh.
			@param mesh is a constant pointer to the mesh to which the Edge belongs
		*/
		Edge_ID(const UInt edge_id, Geometry::Rigid_Mesh * const mesh);

		//! No empty constructor
		Edge_ID () = delete;

		//! Destructor
		virtual ~Edge_ID () {};
		//@}

		//! @name Get Methods
		//@{

		//! Get mesh (const)
		/*!
		 * @return A const pointer to the mesh
		 */
		Geometry::Rigid_Mesh * getMesh () const
			{ return M_mesh; }

		//! Get Edge id (const)
		/*!
		 * @return the edge id contained in the class
		 */
		UInt getEdgeId () const
			{return Edge_Id;}

		//! Get Edge (const)
		/*!
		 * @return a reference to a Edge whose id is contained in the class
		 */
		const Rigid_Mesh::Edge & getEdge () const
			{return M_mesh->getEdgesVector()[Edge_Id];}

		//! Get center (const)
		/*!
		 * @return the centroid of the Edge whose id is contained in the class
		 */
		const Generic_Point getCentroid () const
			{return M_mesh->getEdgesVector()[Edge_Id].getCentroid();}

		//! Get size (const)
		/*!
		 * @return the size of the Edge whose id is contained in the class
		 */
		Real getSize () const
			{return M_mesh->getEdgesVector()[Edge_Id].length();}

		//! Get separated facets (const)
		/*!
		 * @return the id of the separated facets from Edge whose id is contained in the class
		 */
		const std::vector<UInt> & getSeparated () const
			{return M_mesh->getEdgesVector()[Edge_Id].getSeparatedFacetsIds();}
		//@}
	};

	//! Class that represents a border edge
	/*!
		@class Border_Edge
   		This class is derived from the class Edge_ID.
		This class contains the id of an edge which belongs to the border of the domain.
   	*/
	class Border_Edge: public Edge_ID {
	public:
		//! @name Constructor & Destructor
		//@{

		//! Constructor for a Border_Edge given an edge id
		/*!
			@param edge_id is the id of an Edge of a Rigid_Mesh
			@param mesh is a pointer to the mesh to which the edge belongs
		*/
		Border_Edge (const UInt edge_Id, Geometry::Rigid_Mesh * const mesh);

		//! Copy constructor for a Border_Edge given a border_edge belonging to another Rigid_Mesh.
		/*!
		 * @param border_edge reference to a Border_Edge
		 * @param mesh is a pointer to the mesh to which the edge belongs
		 */
		Border_Edge (const Border_Edge & border_edge, Geometry::Rigid_Mesh * const mesh);

		//! Copy constructor for a Border_Edge given a border_edge.
		/*!
		 * @param e reference to a Border_Edge
		 */
		Border_Edge (const Border_Edge & e);

		//! Destructor
		~Border_Edge () = default;
		//@}
	};

	//! Class that represents a fracture edge
	/*!
		@class Fracture_Edge
   		This class is derived from the class Edge_ID.
		This class contains the id of an edge which belongs to at least a fracture.
   	*/
	class Fracture_Edge: public Edge_ID {
	public:
		//! @name Constructor & Destructor
		//@{

		//! Constructor for a Fracture_Edge given an edge id
		/*!
			@param edge_id is the id of an Edge of a Rigid_Mesh
			@param mesh is a pointer to the mesh to which the edge belongs
		*/
		Fracture_Edge (const UInt edge_Id, Geometry::Rigid_Mesh * const mesh);

		//! Copy constructor for a Fracture_Edge given a fracture_edge belonging to another Rigid_Mesh.
		/*!
		 * @param fracture_edge reference to a Fracture_Edge
		 * @param mesh is a pointer to the mesh to which the edge belongs
		 */
		Fracture_Edge (const Fracture_Edge & fracture_edge, Geometry::Rigid_Mesh * const mesh);

		//! Copy constructor for a Fracture_Edge given a fracture_edge.
		/*!
		 * @param e reference to a Fracture_Edge
		 */
		Fracture_Edge (const Fracture_Edge & e);

		//! Destructor
		~Fracture_Edge () = default;
		//@}

		//! @name Get Methods
		//@{

		//! Get fractures ids (const)
		/*!
		 * @return the ids of the represented fractures
		 */
		const std::vector<UInt> & getFractureIds () const
			{return Fracture_Ids;}
		//@}

	protected:
		//! Ids of the represented fractures
		std::vector<UInt> Fracture_Ids;
	};

	//! Class that represents a border fracture edge
	/*!
		@class Border_Fracture_Edge
   		This class is derived from the class Edge_ID.
		This class contains the id of an edge which belongs to the border and belongs to at least a fracture.
   	*/
	class Border_Fracture_Edge: public Edge_ID {
	public:
		//! @name Constructor & Destructor
		//@{

		//! Constructor for a Fracture_Edge given an edge id
		/*!
			@param edge_id is the id of an Edge of a Rigid_Mesh
			@param mesh is a pointer to the mesh to which the edge belongs
		*/
		Border_Fracture_Edge (const UInt edge_Id, Geometry::Rigid_Mesh * const mesh);

		//! Copy constructor for a Fracture_Edge given a fracture_edge belonging to another Rigid_Mesh.
		/*!
		 * @param fracture_edge reference to a Fracture_Edge
		 * @param mesh is a pointer to the mesh to which the edge belongs
		 */
		Border_Fracture_Edge (const Fracture_Edge & fracture_edge, Geometry::Rigid_Mesh * const mesh);

		//! Copy constructor for a Fracture_Edge given a fracture_edge.
		/*!
		 * @param e reference to a Fracture_Edge
		 */
		Border_Fracture_Edge (const Fracture_Edge & e);

		//! Destructor
		~Border_Fracture_Edge () = default;
		//@}

		//! @name Get Methods
		//@{

		//! Get the Border_Edge represented by this edge (const)
		/*!
		 * @return the Border_Edge represented by this edge
		 */
		const Border_Edge & getBorderEdge () const
			{return *M_borderEdge;}

		//! Get the Fracture_Edge represented by this edge (const)
		/*!
		 * @return the Fracture_Edge represented by this edge
		 */
		const Fracture_Edge & getFractureEdge () const
		{return *M_fractureEdge;}
		//@}

		//! @name Set Methods
		//@{

		//! Set the Border_Edge represented by this edge
		/*!
		 * @param the Border_Edge represented by this edge
		 */
		void setBorderEdge(const Border_Edge & borderEdge)
			{ M_borderEdge = &borderEdge;}

		//! Set the Fracture_Edge represented by this edge
		/*!
		 * @param the Fracture_Edge represented by this edge
		 */
		void setFractureEdge(const Fracture_Edge & fractureEdge)
			{ M_fractureEdge = &fractureEdge;}
		//@}

	protected:
		const Border_Edge * M_borderEdge;
		const Fracture_Edge * M_fractureEdge;
	};

	//! Class that handles the several types of facets
	/*!
		@class Facet_ID
   		This class is a container of a Facet-Id.
		This class is a base class.
   	*/
	class Facet_ID{

	protected:
		//! The Id of the referred Facet
		UInt Facet_Id;
		//! Zone code of the facet
		UInt M_zoneCode;
		//! The pointer to the Rigid_Mesh containing the Facet
		Rigid_Mesh * M_mesh;

	public:
		//! @name Constructor & Destructor
		//@{

		//! Constructor for a Facet_ID.
		/*!
			@param facet_id is the id of a Facet of a Rigid_Mesh.
			@param zoneCode the code of the zone
			@param mesh is a constant pointer to the mesh to which the Facet belongs
		*/
		Facet_ID(const UInt facet_id, const UInt zoneCode, Geometry::Rigid_Mesh * const mesh);

		//! No empty constructor
		Facet_ID () = delete;

		//! Destructor
		virtual ~Facet_ID () {};
		//@}

		//! @name Get Methods
		//@{
		
		//! Get mesh (const)
		/*!
		 * @return A const pointer to the mesh
		 */
		Geometry::Rigid_Mesh * getMesh () const
			{ return M_mesh; }
		
		//! Get Facet id (const)
		/*!
		 * @return the facet id contained in the class
		 */
		UInt getFacetId () const
			{return Facet_Id;}
		
		//! Zone Code (const)
		/*!
		 * @return the zone code of the represented fractures
		 */
		UInt getZoneCode () const
			{return M_zoneCode;}

		//! Get Facet (const)
		/*!
		 * @return a reference to a Facet whose id is contained in the class
		 */
		const Rigid_Mesh::Facet & getFacet () const
			{return M_mesh->getFacetsVector()[Facet_Id];}
		
		//! Get unsigned normal (const)
		/*!
		 * @return the unsigned normal vector to the Facet whose id is contained in the class
		 */
		const Generic_Vector getUNormal () const
			{return M_mesh->getFacetsVector()[Facet_Id].getUnsignedNormal();}
		
		//! Get center (const)
		/*!
		 * @return the centroid of the Facet whose id is contained in the class
		 */
		const Generic_Point getCentroid () const
			{return M_mesh->getFacetsVector()[Facet_Id].getCentroid();}
		
		//! Get size (const)
		/*!
		 * @return the size of the Facet whose id is contained in the class
		 */
		Real getSize () const
			{return M_mesh->getFacetsVector()[Facet_Id].area();}
		
		//! Get separated cells (const)
		/*!
		 * @return the id of the separated cell from Facet whose id is contained in the class
		 */
		const std::vector<UInt> & getSeparated () const
			{return M_mesh->getFacetsVector()[Facet_Id].getSeparatedCellsIds();}
		//@}
	};

public:
	
	//! Class that represents an interior facet (no fracture facet)
	/*!
		@class Regular_Facet
   		This class is derived from the class Facet_ID.
		This class contains the id of a Facet which is neither a border facet nor a fracture facet.
   	*/
	class Regular_Facet : public Facet_ID {
	public:
		//! @name Constructor & Destructor
		//@{

		//! Constructor for a Regular_Facet given an id
		/*!
			@param facet_id is the id of a Facet of a Rigid_Mesh
			@param zoneCode the code of the zone
			@param mesh is a pointer to the mesh to which the Facet belongs
		*/
		Regular_Facet (const UInt facet_Id, const UInt zoneCode, Geometry::Rigid_Mesh * const mesh);

		//! Copy constructor for a Regular_Facet given a regular_facet of another Rigid_Mesh.
		/*!
			@param regular_facet is a Regular_Facet of a different Rigid_Mesh
			@param mesh is a pointer to the mesh to which the Facet belongs
		*/
		Regular_Facet (const Regular_Facet & regular_facet, Geometry::Rigid_Mesh * const mesh);

		//! Copy constructor for a Regular_Facet given a regular_facet
		/*!
		 * @param f reference to a Regular_Facet
		 */
		Regular_Facet (const Regular_Facet & f);

		//! Destructor
		~Regular_Facet () = default;
		//@}
	};

	//! Class that represents a border facet
	/*!
		@class Border_Facet
   		This class is derived from the class Facet_ID.
		This class contains the id of a facet which belongs to the border of the domain.
   	*/
	class Border_Facet: public Facet_ID {
	protected:
		//! Identifier of the border
		UInt Border_Id;
	public:
		//! @name Constructor & Destructor
		//@{

		//! Constructor for a Border_Facet given a facet id and the relative border id
		/*!
			@param facet_id is the id of a Facet of a Rigid_Mesh
			@param zoneCode the code of the zone
			@param border_id is the id of the border to which the facet belongs.
			@param mesh is a pointer to the mesh to which the facet belongs
		*/
		Border_Facet (const UInt facet_Id, const UInt zoneCode, const UInt border_Id, Geometry::Rigid_Mesh * const mesh);

		//! Copy constructor for a Border_Facet given a border_facet belonging to another Rigid_Mesh.
		/*!
		 * @param border_facet reference to a Border_Facet
		 * @param mesh is a pointer to the mesh to which the facet belongs
		 */
		Border_Facet (const Border_Facet & border_facet, Geometry::Rigid_Mesh * const mesh);

		//! Copy constructor for a Border_Facet given a border_facet.
		/*!
		 * @param f reference to a Border_Facet
		 */
		Border_Facet (const Border_Facet & f);

		//! Destructor
		~Border_Facet () = default;
		//@}

		//! @name Get Methods
		//@{
		//! Get border id (const)
		/*!
		 * @return the id of the border of the contained Facet
		 */
		UInt getBorderId() const
			{return Border_Id;}
		//@}
	};

	//! Class that represents a fracture facet
	/*!
		@class Fracture_Facet
   		This class is derived from the class Facet_ID.
		This class contains the id of a facet which belongs to a fracture.
   	*/
	class Fracture_Facet: public Facet_ID {
	public:
		//! @name Constructor & Destructor
		//@{

		//! Constructor for a Fracture_Facet given a facet id and the id of the represented fractures.
		/*!
			@param facet_Ids is a triplet: 1. id of a Facet of a Rigid_Mesh 2. id as Fracture_Facet 3. number of cells in Rigid_Mesh 4. Zone Code
			@param fracture_Ids is a vector with the ids of the fractures represented by the Facet.
			@param mesh is a pointer to the mesh to which the facet belongs
		*/
		Fracture_Facet(const std::tuple<UInt,UInt,UInt,UInt> facet_Ids, const std::set<UInt> & fracture_Ids, Geometry::Rigid_Mesh * const mesh);

		//! Copy constructor for a Fracture_Facet given a fracture_facet belonging to another Rigid_Mesh.
		/*!
		 * @param fracture_facet reference to a Fracture_Facet
		 * @param mesh is a pointer to the mesh to which the facet belongs
		 */
		Fracture_Facet(const Fracture_Facet & fracture_facet, Geometry::Rigid_Mesh * const mesh);

		//! Copy constructor for a Fracture_Facet
		/*!
		 * @param fracture_facet reference to a Fracture_Facet
		 */
		Fracture_Facet(const Fracture_Facet&);

		//! Destructor
		~Fracture_Facet () = default;
		//@}

		//! @name Get Methods
		//@{

		//! Get the id as Fracture_Facet (const)
		/*!
		 * @return the id as Fracture_Facet
		 */
		UInt getId () const
			{return M_Id;}

		//! Get the number of cells (const)
		/*!
		 * @return the number of cells in the Rigid_Mesh
		 */
		UInt getCellsnumber () const
			{return Cells_number;}

		//! Get the id as Cell (const)
		/*!
		 * @return the id as Cell = Number of cells + Id as Fracture_Facet
		 */
		UInt getIdasCell () const
			{return (M_Id+Cells_number);}

		//! Get fractures ids (const)
		/*!
		 * @return the ids of the represented fractures
		 */
		const std::vector<UInt> & getFractureIds () const
			{return Fracture_Ids;}

		//! Get fracture neighbors (const)
		/*!
		 * @return a map that associate a juncture the ids of the neighboring Fracture_Facet
		 */
		const std::map<Fracture_Juncture, std::vector<UInt> > & getFractureNeighbors () const
			{return Fracture_Neighbors;}

		//! Get fracture tips (const)
		/*!
		 * @return a set of tips of the current fracture
		 */
		const std::set<Fracture_Tip> & getFractureTips () const
			{return Fracture_Tips;}
		//@}

		friend class Rigid_Mesh;

	private:
		//! Number of cells in the Rigid_Mesh
		UInt Cells_number;
		//! Id as Fracture_Facet
		UInt M_Id;
		//! Ids of the represented fractures
		std::vector<UInt> Fracture_Ids;
		//! Map that associates a juncture the ids of the neighboring Fracture_Facet
		std::map<Fracture_Juncture, std::vector<UInt> > Fracture_Neighbors;
		//! Set of tips of the current fracture
		std::set<Fracture_Tip> Fracture_Tips;
	};

public:
	//! @name Constructor & Destructor
	//@{

	//! Constructor for a Rigid_Mesh given a Generic_Mesh and a Geometry::PropertiesMap
	/*!
		@param generic_mesh is a reference to a Generic_Mesh from which the Rigid_Mesh is constructed
		@param prop reference to a Geometry::PropertiesMap
		@param renumber If False the facets and cells are not renumbered (Default = false)
	*/
	Rigid_Mesh (Generic_Mesh & generic_mesh, const PropertiesMap & prop, const bool renumber = false);

	//! Copy constructor
	/*!
	 * @param mesh reference to a Rigid_Mesh
	 */
	Rigid_Mesh (const Rigid_Mesh & mesh);

	//! No empty constructor
	Rigid_Mesh () = delete;

	//! Destructor
	~Rigid_Mesh () = default;

	//@}

	//! @name Get Methods
	//@{

	//! Get nodes vector (const)
	/*!
	 * @return a reference to the vector that contains the nodes of the mesh
	 */
	const std::vector<Point> & getNodesVector () const
		{ return M_nodes; }

	//! Get edges vector (const)
	/*!
	 * @return a reference to the vector that contains the edges of the mesh
	 */
	const std::vector<Edge> & getEdgesVector () const
		{ return M_edges; }

	//! Get facets vector (const)
	/*!
	 * @return a reference to the vector that contains the facets of the mesh
	 */
	const std::vector<Facet> & getFacetsVector () const
		{ return M_facets; }

	//! Get cells vector (const)
	/*!
	 * @return a reference to the vector that contains the cells of the mesh
	 */
	const std::vector<Cell> & getCellsVector () const
		{ return M_cells; }

	//! Get border edges ids vector (const)
	/*!
	 * @return a reference to the vector of the Border_Edge
	 */
	const std::vector<Border_Edge> & getBorderEdgesIdsVector () const
		{ return M_borderEdges; }

	//! Get internal facets ids vector (const)
	/*!
	 * @return a reference to the vector of the Regular_Facet
	 */
	const std::vector<Regular_Facet> & getInternalFacetsIdsVector () const
		{ return M_internalFacets; }

	//! Get border facets ids vector (const)
	/*!
	 * @return a reference to the vector of the Border_Facet
	 */
	const std::vector<Border_Facet> & getBorderFacetsIdsVector () const
		{ return M_borderFacets; }

	//! Get fracture facets ids vector (const)
	/*!
	 * @return a reference to the vector of the Fracture_Facet
	 */
	const std::vector<Fracture_Facet> & getFractureFacetsIdsVector () const
		{ return M_fractureFacets; }

	//! Get size of biggest facet (const)
	/*!
	 * @return the size of the biggest facet
	 */
	Real getMaxFacetSize () const
	{ 
		auto comp=[](Facet f1, Facet f2){return f1.area()<f2.area();};
		return std::max_element(M_facets.begin(),M_facets.end(), comp)->area();
	}

	//! Get size of smallest facet (const)
	/*!
	 * @return the size of the smallest facet
	 */
	Real getMinFacetSize () const
	{ 
		auto comp=[](Facet f1, Facet f2){return f1.area()<f2.area();};
		return std::min_element(M_facets.begin(),M_facets.end(), comp)->area();
	}

	//! Get average size of facet (const)
	/*!
	 * @return the average size of the facets
	 */
	Real getAveFacetSize () const
	{ 	
		Real ave = 0;
		auto sum = [&ave](Facet f1){ave+=f1.area();};
		for(auto it: M_facets) sum(it);
		return ave/M_facets.size();
	}

	//! Get the properties map (const)
	/*!
	 * @return the properties of the matrix and fractures
	 */
	const PropertiesMap & getPropertiesMap () const
		{ return M_properties; }

	//! Are facets and cells renumbered? (const)
	/*!
	 * @return True if the facets and cells are renumbered
	 */
	bool isRenumbered () const
		{ return M_renumber; }
	//@}

	//! Test if a facet has as neighbor a specific cell
	/*!
	 * @param facet_Id is the Facet Id
	 * @param idNeighbour is the Cell Id
	 * @return True if the facet has idNeighbour as neighbor through Facet facet_Id
	 */
	bool hasNeighborsThroughFacet (const UInt & facet_Id, const UInt & idNeighbor) const;

	//! @name Methods
	//@{

	//! Display general information about the Rigid_Mesh
	/*!
	 * @param out specify the output format (std::cout by default)
	 */	
	void showMe ( std::ostream & out=std::cout ) const;

	//! Converts ids of nodes to nodes
	/*!
	 * @param pointsId a reference to a vector of ids of nodes
	 * @return a vector with the corresponding points
	*/
	const std::vector<Generic_Point> IdToPoints (const std::vector<UInt> & pointsId);

	//@}

protected:
	//! @name Protected Methods
	//@{

	//! It is called by the constructor
	/*!
	 * @param generic_mesh the generic_mesh
	 */
	void M_constructor (Generic_Mesh & generic_mesh);

	//! Prints the components of a point
	/*!
	 * @param generic_point is the point to print
	 * @param out specify the output format (std::cout by default)
	 */	
	void showPoint (const Point & generic_point, std::ostream & out=std::cout ) const;

	//! Builds the vector of cells. It is called by constructor
	/*!
	 * @param generic_mesh the generic_mesh
	 * @param old_to_new_map is a map that binds the id of a cell in the original Generic_Mesh to the id in the Rigid_Mesh
	 */	
	void CellsVectorBuilder ( Generic_Mesh & generic_mesh, std::map<UInt, UInt> & old_to_new_map );

	//! Builds the vector of facets. It is called by constructor
	/*!
	 * @param generic_mesh the generic_mesh
	 * @param old_to_new_mapCells is a map that binds the id of a cell in the original Generic_Mesh to the id in the Rigid_Mesh
	 * @param old_to_new_mapFacets is a map that binds the id of a facet in the original Generic_Mesh to the id in the Rigid_Mesh
	 */	
	void FacetsVectorsBuilder ( Generic_Mesh & generic_mesh, const std::map<UInt, UInt> & old_to_new_mapCells, std::map<UInt, UInt> & old_to_new_mapFacets);

	//! Builds the vector of edges. It is called by constructor
	void EdgesVectorBuilder();

	//! Builds the vector of facets. It is called by constructor
	/*!
	 * @param old_to_new_mapCells a map that binds the id of a cell in the original Generic_Mesh to the id in the Rigid_Mesh
	 * @param old_to_new_mapFacets a map that binds the id of a facet in the original Generic_Mesh to the id in the Rigid_Mesh
	 * @param generic_mesh reference to the Generic_Mesh
	*/
	void M_facetsVectorsBuilder(const std::map<UInt, UInt> & old_to_new_mapCells, std::map<UInt, UInt> & old_to_new_mapFacets, Generic_Mesh & generic_mesh);

	//! Converts the neighbors ids of cells from old one to new one. It is called by constructor
	/*!
	 * @param old_to_new_map a map that binds the id of a cell in the original Generic_Mesh to the id in the Rigid_Mesh
	*/
	void AdjustCellNeighbors(const std::map<UInt, UInt> & old_to_new_map);

	//! Converts the facets id of the cells from old one to new one. It is called by constructor
	/*!
	 * @param old_to_new_mapFacets a map that binds the id of a facet in the original Generic_Mesh to the id in the Rigid_Mesh
	 */
	void AdjustCellFacets(const std::map<UInt, UInt> & old_to_new_mapFacets);

protected:
	
	//! Vector of Nodes
	std::vector<Point> M_nodes;
	//! Vector of Edges
	std::vector<Edge> M_edges;
	//! Vector of Facets
	std::vector<Facet> M_facets;
	//! Vector of Cells
	std::vector<Cell> M_cells;

	//! Vector of Border_Edge
	std::vector<Border_Edge> M_borderEdges;
	//! Vector of Fracture_Edge
	std::vector<Fracture_Edge> M_fractureEdges;
	//! Vector of Border_Fracture_Edge
	std::vector<Border_Fracture_Edge> M_borderFractureEdges;

	//! Vector of Regular_Facet
	std::vector<Regular_Facet> M_internalFacets;
	//! Vector of Border_Facet
	std::vector<Border_Facet> M_borderFacets;
	//! Vector of Fracture_Facet
	std::vector<Fracture_Facet> M_fractureFacets;

	//! Contains the properties of the matrix and fractures
	const PropertiesMap & M_properties;

	//! True if the ids are renumbered
	bool M_renumber;
};

}//namespace Geometry

#endif
