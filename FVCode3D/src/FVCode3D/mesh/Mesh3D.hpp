 /*!
 * @file Mesh3D.hpp
 * @brief Classes that implements polyhedrical unstructured meshes.
 */

#ifndef MESH3D_HPP_
#define MESH3D_HPP_

#include <FVCode3D/core/TypeDefinition.hpp>
#include <FVCode3D/geometry/CoordinateSystem.hpp>
#include <FVCode3D/geometry/Operations.hpp>
#include <FVCode3D/fracture/Fracture3D.hpp>
#include <FVCode3D/fracture/FractureNetwork3D.hpp>
#include <FVCode3D/mesh/TetGenWrapper.hpp>

#include <set>
#include <vector>
#include <map>
#include <cmath>

namespace FVCode3D
{

class FractureNetwork3D;
class Fracture;

//! Class that implements a 3D mesh.
/*!
	@class Mesh3D
	This class implements the concept of a three dimensional polyhedrical unstructured mesh.
	It is composed of points, facets and cells.
 */
class Mesh3D
{

public:

	class Facet3D;

public:

	//! Class that implements a facet
	/*!
		@class Facet3D
		This class implements the concept of a facet of a mesh.
		It is composed of a vector of ids that represent the points of the mesh.
	*/
	class Facet3D{

	public:
		//! @name Constructor & Destructor
		//@{

		//! Empty constructor
		Facet3D();

		//! Copy constructor
		/*!
		 * @param facet a reference to a Mesh3D::Facet3D
		 */
		Facet3D(const Facet3D & facet);

		//! Constructor for the facet of a specified mesh moving from the vertices.
		/*!
		 * The normal is computed in clockwise rotation.
		 * @param mesh pointer to a Mesh3D
		 * @param vertices vector of identifiers of the vertices
		 * @param borderId id of the border. If the facet is interior, then borderId should be 0
		 * @pre vertices[i] is connected to vertices[i-i] and vertices[i+i]
		 */
		Facet3D(Mesh3D * const mesh, const std::vector<UInt> & vertices, const UInt zone, const UInt borderId);

		//@}

		//! @name Get Methods
		//@{

		//! Get mesh (const)
		/*!
		 * @return the pointer to the mesh
		 */
		const Mesh3D * getMesh() const
			{ return M_mesh; }

		//! Get the i-th vertex id of the current facet (const)
		/*!
		 * @param i i-th vertex of the current facet
		 * @return the id of i-th vertex
		 */
		UInt getVertexId(const UInt i) const
			{ return M_vertexIds[i]; }

		//! Get the vertex vector of the current facet (const)
		/*!
		 * @return a constant reference to the vertex vector of the current facet
		 */
		const std::vector<UInt> & getVerticesVector() const
			{ return M_vertexIds; }

		//! Get the i-th vertex (const)
		/*!
		 * @return a constant reference to the i-th vertex
		 */
		const Point3D & getVertex(const UInt i) const
			{ return M_mesh->getNodesVector()[M_vertexIds[i]]; }

		//! Get the vertex vector(const)
		/*!
		 * @return a constant reference to the vertex vector
		 */
		const std::vector<Point3D> & getVertices() const
			{ return M_mesh->getNodesVector(); }

		//! Get the set of the separated cells ids (const)
		/*!
		 * @return a constant reference to the set that contains the id of the separated cells
		 */
		const std::set<UInt> & getSeparatedCells() const
			{ return M_separatedCells; }

		//! Get the set of the separated cells ids
		/*!
		 * @return a reference to the set that contains the id of the separated cells
		 */
		std::set<UInt> & getSeparatedCells()
			{ return M_separatedCells; }

		//! Get the set of the represented fractures (const)
		/*!
		 * @return a constant reference to the set that contains the id of the represented fractures
		 */
		const std::set<UInt> & getRepresentedFractureIds() const
			{ return M_representedFractureIds; }

		//! Get the set of the represented fractures
		/*!
		 * @return a reference to the set that contains the id of the represented fractures
		 */
		std::set<UInt> & getRepresentedFractureIds()
			{ return M_representedFractureIds; }

		//! Get the number of points of the facet
		/*!
		 * @return the number of points of the facet
		 */
		UInt getNumberOfVertices() const
			{ return M_vertexIds.size(); }

		//! Get the facet centroid
		/*!
		 * @return a reference to the facet centroid
		 */
		const Point3D & getCentroid() const
			{ return M_centroid; }

		//! Get the borderID
		/*!
		 * @return the borderID. If zero, then the facet is an interior facet
		 */
		UInt getBorderId() const
			{ return M_borderId; }

		//! Get the zone code
		/*!
		 * @return the zone code of the facet
		 */
		UInt getZoneCode() const
			{ return M_zone; }

		//@}

		//! @name Set Methods
		//@{

		//! Set the mesh pointer
		/*!
		 * @param mesh the new value for the mesh pointer
		 */
		void setMesh(Mesh3D * const mesh)
			{ M_mesh = mesh; }

		//! Set the ids of the two vertexes
		/*!
		 * @param vertices the ids of a vertices
		 */
		void setVertices(const std::vector<UInt> & vertices)
			{ M_vertexIds = vertices; }

		//! Set the borderID
		/*!
		 * @param The borderID. If zero, then the facet is an interior facet
		 */
		void setBorderID(const UInt borderId)
			{ M_borderId = borderId; }

		//! Set the zone code
		/*!
		 * @param zone the zone code of the facet
		 */
		void setZoneCode(const UInt zone)
			{ M_zone = zone; }

		//@}

		//! @name Methods
		//@{

		//! Test if the facet represents a fracture
		/*!
		 * @return True if the facet represents a fracture
		 */
		bool isFracture() const
			{ return !M_representedFractureIds.empty(); }

		//! Test if the facet is a border facet
		/*!
		 * @return True if the facet is a border facet (i.e. border id != 0)
		 */
		bool isBorderFacet() const
			{ return M_borderId; }

		//! Compute the normal of the facet
		/*!
		 * @return The normal of the facet
		 */
		Point3D computeNormal() const;

		//! Compute the facet centroid and stores it in the M_centroid variable
		void computeCentroid();

		//! Compute the area of the facet
		/*!
		 * @return the area of the facet
		 */
		Real area() const;

		//! Display general information about the content of the class
		/*!
		 * @param out specify the output format (std::cout by default)
		 */
		void showMe(std::ostream & out=std::cout) const;

		//@}

	private:
		//! The pointer to the mesh containing this facet
		const Mesh3D * M_mesh;
		//! The ids of the vertices of the facet
		std::vector<UInt> M_vertexIds;
		//! The set containing the ids of the cells separated by this facet
		std::set<UInt> M_separatedCells;
		//! The set containing the ids of the represented fractures
		std::set<UInt> M_representedFractureIds;
		//! The centroid of the facet
		Point3D M_centroid;
		//! Border Id. If zero, then the facet is an interior facet
		UInt M_borderId;
		//! Zone code
		UInt M_zone;
	};

	//! Class that implements a cell
	/*!
		@class Cell3D
		This class implements the concept of a Cell of a mesh.
		It is composed of a vector of ids that represent the facets of the mesh.
	 */
	class Cell3D{

	public:
		//! @name Constructor & Destructor
		//@{

		//! Empty constructor
		Cell3D();

		//! Copy constructor
		/*!
		 * @param cell a reference to a Mesh3D::Cell3D
		 */
		Cell3D(const Cell3D & cell);

		//! Constructor for a cell in a specified mesh given the facets
		/*!
		 * @param mesh pointer to the mesh
		 * @param facets the vector that contains the facets ids
		 * @param zone code of the zone
		 */
		Cell3D( const Mesh3D * mesh,
				const std::vector<UInt> & facets,
				const UInt zone );

		//@}

		//! @name Get Methods
		//@{

		//! Get mesh (const)
		/*!
		 * @return the pointer to the mesh
		 */
		const Mesh3D * getMesh() const
			{ return M_mesh; }

		//! Get the vector that contains the vertexes ids (const)
		/*!
		 * @return a constant reference to the vector of the vertexes ids
		 */
		const std::vector<UInt> & getVerticesVector() const
			{ return M_vertexIds; }

		//! Get the vector that contains the vertexes ids
		/*!
		 * @return a reference to the vector of the vertexes ids
		 */
		std::vector<UInt> & getVerticesVector()
			{ return M_vertexIds; }

		//! Get the set that contains the facets ids (const)
		/*!
		 * @return a constant reference to the set of the facets ids
		 */
		const std::set<UInt> & getFacetsSet() const
			{ return M_facetIds; }

		//! Get the set that contains the facets ids
		/*!
		 * @return a reference to the set of the facets ids
		 */
		std::set<UInt> & getFacetsSet()
			{ return M_facetIds; }

		//! Get the set that contains the ids of the neighbor cells (const)
		/*!
		 * @return a constant reference to the set of the ids of the neighbor cells
		 */
		const std::set<UInt> & getNeighborsSet() const
			{ return M_neighborIds; }

		//! Get the set that contains the ids of the neighbor cells
		/*!
		 * @return a reference to the set of the ids of the neighbor cells
		 */
		std::set<UInt> & getNeighborsSet()
			{ return M_neighborIds; }

		//! Get the cell centroid
		/*!
		 * @return a reference to the cell centroid
		 */
		const Point3D & getCentroid() const
			{ return M_centroid; }

		//! Return the volume (const)
		/*!
		 * @return the volume of the cell
		 */
		Real volume() const
			{ return M_volume; }

		//! Get the zone code
		/*!
		 * @return the zone code of the cell
		 */
		UInt getZoneCode() const
			{ return M_zone; }

		//@}

		//! @name Set Methods
		//@{

		//! Set the zone code
		/*!
		 * @param zone the zone code of the cell
		 */
		void setZoneCode(const UInt zone)
			{ M_zone = zone; }

		//@}

		//! @name Methods
		//@{

		//! The number of vertexes of the cell
		/*!
		 * @return the number of vertexes
		 */
		UInt verticesNumber() const
			{ return M_vertexIds.size(); }

		//! The number of facets of the cell
		/*!
		 * @return the number of facets
		 */
		UInt facetsNumber() const
			{ return M_facetIds.size(); }

		//! The number of neighbors of the cell
		/*!
		 * @return the number of neighbors
		 */
		UInt neighborsNumber() const
			{ return M_neighborIds.size(); }

		//! The normalized vector normal to the facet defined by the facet id
		/*!
		 * @param facetId the facet id
		 * @return The normalized vector normal to the facet
		 */
		Point3D outerNormalToFacet( const UInt & facetId) const;

		//! Test if a cell has a neighbor cell through the facet defined by the facet id
		/*!
		 * @param facetId The facet id
		 * @param idNeighbor The variable in which the id of the neighboring cell is stored
		 * @return TRUE  -> the neighboring cell exist
		 *		   FALSE -> the neighboring cell does not exist
		 */
		bool hasNeighborsThroughFacet( const UInt & facetId, UInt & idNeighbor) const;

		//! Computes the volume and the centroid of the cell
		/*!
		 * @return the volume and the ccentroid of the cell
		 */
		void computeVolumeAndCentroid();

		//! Display general information about the content of the class
		/*!
		 * @param out specify the output format (std::cout by default)
		 */
		void showMe(std::ostream & out=std::cout) const;

		//@}

	private:

		//! The pointer to the mesh containing this cell
		const Mesh3D * M_mesh;
		//! The vector containing the vertexes ids. The order doesn't matter.
		std::vector<UInt> M_vertexIds;
		//! The vector containing the facets ids
		std::set<UInt> M_facetIds;
		//! The set storing the ids of the neighbor cells
		std::set<UInt> M_neighborIds;
		//! The centroid of the cell
		Point3D M_centroid;
		//! The volume of the cell
		Real M_volume;
		//! Zone code
		Real M_zone;
	};

public:

	//! @name Constructor & Destructor
	//@{

	//! Empty constructor
	Mesh3D();

	//! Copy constructor
	/*!
	 * @param mesh reference to a Mesh3D
	 */
	Mesh3D(const Mesh3D & mesh);

	//! Destructor
	virtual ~Mesh3D() = default;

	//@}

	//! @name Get Methods
	//@{

	//! Get the vector storing the nodes of the mesh (const)
	/*!
	 * @return a constant reference to the vector storing the nodes of the mesh
	 */
	const std::vector<Point3D> & getNodesVector() const
		{ return M_nodes; }

	//! Get the vector storing the nodes of the mesh
	/*!
	 * @return a reference to the vector storing the nodes of the mesh
	 */
	std::vector<Point3D> & getNodesVector()
		{ return M_nodes; }

	//! Get the set storing the facets of the mesh (const)
	/*!
	 * @return a constant reference to the set storing the facets of the mesh
	 */
	const std::map<UInt,Facet3D> & getFacetsMap() const
		{ return M_facets; }

	//! Get the map storing the facets of the mesh
	/*!
	 * @return a reference to the map storing the facets of the mesh
	 */
	std::map<UInt,Facet3D> & getFacetsMap()
		{ return M_facets; }

	//! Get the map storing the cells of the mesh (const)
	/*!
	 * @return a constant reference to the map storing the cells of the mesh
	 */
	const std::map<UInt,Cell3D> & getCellsMap() const
		{ return M_cells; }

	//! Get the map storing the cells of the mesh
	/*!
	 * @return a reference to the map storing the cells of the mesh
	 */
	std::map<UInt,Cell3D> & getCellsMap()
		{ return M_cells; }

	//! Get the fracture network of the mesh (const)
	/*!
	 * @return a constant reference to the fracture network of the mesh
	 */
	const FractureNetwork3D & getFN() const
		{ return M_fn; }

	//! Get the fracture network of the mesh
	/*!
	 * @return a constant reference to the fracture network of the mesh
	 */
	FractureNetwork3D & getFN()
		{ return M_fn; }

	//@}

	//! @name Methods
	//@{

	//! Get the number of nodes
	UInt getNumberOfNodes() const
		{ return M_nodes.size(); }

	//! Get the number of facets
	UInt getNumberOfFacets() const
		{ return M_facets.size(); }

	//! Get the number of cells
	UInt getNumberOfCells() const
		{ return M_cells.size(); }

	//! Add the fracture network
	void addFractureNetwork( FractureNetwork3D & fn )
		{ M_fn.getNetwork() = fn.getNetwork(); }

	//! Update the Facet3D with the id of the fracture
	/*!
	 * @pre addFractureNetwork()
	 */
	void updateFacetsWithFractures();

	//! Update the Facet3D with the id of the separated cells
	void updateFacetsWithCells();

	//! Update the Cell3D with the neighboring cells
	/*!
	 * @pre updateFacetsWithCells()
	 */
	void updateCellsWithNeighbors();

	//! Fill a map that associates to node ids of each facet the facet id
	void buildNodesToFacetMap();

	//! Find the facet id given the list of the nodes
	/*!
	 * @param nodes vector of nodes that define the facet
	 * @return the id of the facet
	 */ 
	UInt getFacetFromNodes(std::vector<UInt> & nodes);
	
	//! Export the mesh in vtu format
	/*!
	 * (Use paraview to open the vtu file)
	 * @param filename The name of the vtu file created by this method
	 * @return TRUE  -> operation ended correctly
	 *		   FALSE -> an error occurred
	 */
	bool exportVTU(const std::string & filename) const;

	//! Export some cells of the mesh in vtu format
	/*!
	 * (Use paraview to open the vtu file)
	 * @param filename The name of the vtu file created by this method
	 * @return TRUE  -> operation ended correctly
	 *		   FALSE -> an error occurred
	 */
	bool exportCellsVTU(const std::string & filename, const std::vector<UInt> & idCells) const;

	//! Export in a single vtk file all the facets representing each fracture in the network
	/*!
	 * (Use paraview to open the vtk file)
	 * @param filename The name of the vtk file created by this method
	 * @return TRUE  -> operation ended correctly
	 *		   FALSE -> an error occurred
	 */
	bool exportFractureNetworkVTK(const std::string & filename) const;

	//! Export in a vtk file the facets representing a fracture in the network
	/*!
	 * (Use paraview to open the vtk file)
	 * @param filename The name of the vtk file created by this method
	 * @param f the fracture id
	 * @return TRUE  -> operation ended correctly
	 *		   FALSE -> an error occurred
	 */
	bool exportFractureVTK(const std::string & filename, const UInt & f) const;

    //! Clear the mesh
    void clear();

	//@}

protected:
	//! The fracture network
	FractureNetwork3D M_fn;
	//! The vector storing all the nodes of the mesh
	std::vector<Point3D> M_nodes;
	//! The set storing all the edges of the mesh
	std::map<UInt,Facet3D> M_facets;
	//! The map storing all the cells of the mesh
	std::map<UInt,Cell3D> M_cells;
	//! The map that associates to a vector of node ids a facet id
	std::map<std::vector<UInt>, UInt> M_nodesToFacet;

};

//! Overloading of operator< between two Facet3D objects (needed by the set M_facets)
/*!
 * @param f1 The first facet to be tested
 * @param f2 The second facet to be tested
 * @return TRUE  -> f1 < f2
 *		   FALSE -> f1 >= f2
 */
bool operator<(const Mesh3D::Facet3D & f1, const Mesh3D::Facet3D & f2);

} // namespace FVCode3D

#endif /* MESH3D_HPP_ */
