 /*!
 *	@file Mesh2D.hpp
 *	@brief Class for unstructured mesh in 2D space.
 *
 *	@author Luca Turconi <luca.turconi@moxoff.com>
 *
 */ 

#ifndef MESH2D_HPP_
#define MESH2D_HPP_
 
#include "TypeDefinition.hpp"
#include "Fracture2D.hpp"
#include "Triangulation2D.hpp"

#include<set>
#include<vector>
#include<map>


namespace Geometry{

	/*!
		@class Mesh2D
    	This class implements the concept of mesh in 2D space.
    */
class Mesh2D{
public:

		/*!
			@class Edge2D
    		This class implements the concept of Edge of a Mesh2D.
    	*/
	class Edge2D{
	public:
	
		//! @name Constructor & Destructor
		//@{
	
		//! Empty constructor
		Edge2D();
	
		//! Constructor for an empty edge of a specified mesh
		Edge2D(Geometry::Mesh2D * const mesh);
	
		//! Copy constructor
		Edge2D(const Edge2D & edge);
	
		//! Constructor for the edge (p,q) of a specified mesh
		Edge2D(Geometry::Mesh2D * const mesh, const UInt & p, const UInt & q);
	
		//@}
	
		//! @name Get Methods
		//@{
		
		//! Get mesh (const)
		/*!
		 * @return The pointer to the related mesh
		 */
		Geometry::Mesh2D * const getMesh() const
			{ return M_mesh; }
	
		//! Get the first vertex id (const)
		/*!
		 * @return A constant reference to the first vertex id
		 */
		const UInt & getP1() const
			{ return M_idP1; }

		//! Get the second vertex id (const)
		/*!
		 * @return A constant reference to the second vertex id
		 */
		const UInt & getP2() const
			{ return M_idP2; }
		
		//! Get the first vertex (const)
		/*!
		 * @return A constant reference to the first vertex
		 */
		const Geometry::Point2D & getPointP1() const;
		
		//! Get the second vertex (const)
		/*!
		 * @return A constant reference to the second vertex
		 */
		const Geometry::Point2D & getPointP2() const;
	
		//! Get the segment representing the edge
		/*!
		 * @return The segment between the first and the second vertex
		 */
		Geometry::Segment2D segment() const;

		//! Get the set of the separated cells ids (const)
		/*!
		 * @return A constant reference to the set containing the id of the separated cells
		 */
		const std::set<UInt> & getSeparatedCells() const
			{ return M_separatedCells; }

		//! Get the set of the separated cells ids
		/*!
		 * @return A reference to the set containing the id of the separated cells
		 */
		std::set<UInt> & getSeparatedCells()
			{ return M_separatedCells; }
		
		//! Get the set of the represented fractures (const)
		/*!
		 * @return A constant reference to the set containing the id of the represented fractures
		 */
		const std::set<UInt> & getRepresentedFractureIds() const
			{ return M_representedFractureIds; }
			
		//! Get the set of the represented fractures
		/*!
		 * @return A reference to the set containing the id of the represented fractures
		 */			
		std::set<UInt> & getRepresentedFractureIds()
			{ return M_representedFractureIds; }
		
		//@}
	
		//! @name Set Methods
		//@{
		
		//! Set the mesh pointer
		/*!
		 * @param mesh The new value for the mesh pointer
		 */	
		void setMesh(Geometry::Mesh2D * const mesh)
			{ M_mesh = mesh; }

		//! Set the ids of the two vertexes
		/*!
		 * @param p The id of a vertex
		 * @param q The id of another vertex
		 */		
		void setP(const UInt & p, const UInt & q);
	
		//@}
	
		//! @name Methods
		//@{
		
		//! Test if the edge represent some fracture
		/*!
		 * @return True if the edge represent some fracture
		 */
		bool isFracture() const
			{ return !M_representedFractureIds.empty(); }
		
		//! Test if the edge is a border edge
		/*!
		 * @return True if the edge is a border edge
		 */
		bool isBorderEdge() const
			{ return M_separatedCells.size()==1 ? 1 : 0 ; }
		
		//! Compute the length of the edge
		/*!
		 * @return The length of the edge
		 */		
		Real length() const;
		
		//! Compute the average between the length of the edges of the equilateral triangles with the same area of the separated cells
		/*!
		 * @return The average of the equilateral triangles edges
		 */		
		Real equilateralLength() const;
	
		//! Display general information about the content of the class
		/*!
		 * @param out Specify the output format (std::cout by default)
		 */	
		void showMe(std::ostream & out=std::cout) const;
	
		//@}
	
	private:
		//! The pointer to the mesh containing this edge
		Geometry::Mesh2D * M_mesh;
		//! The id of the first vertex of the edge (M_idP1<M_idP2)
		UInt M_idP1;
		//! The id of the second vertex of the edge (M_idP2>M_idP1)
		UInt M_idP2;
		//! The set containing the ids of the cells separated by this edge
		std::set<UInt> M_separatedCells;
		//! The set containing the ids of the represented fractures
		std::set<UInt> M_representedFractureIds;
		
	};

public:
		/*!
			@class Cell2D
    		This class implements the concept of Cell of a Mesh2D.
    	*/
	class Cell2D{
	public:
		//! @name Constructor & Destructor
		//@{
	
		//! Empty constructor
		Cell2D();
		
		//! Constructor for an empty cell in a specified mesh
		/*!
		 * @param mesh The mesh
		 */
		Cell2D(Geometry::Mesh2D * const mesh);
	
		//! Copy constructor
		Cell2D(const Cell2D & cell);
	
		//! Constructor for a cell in a specified mesh given the vertexes and the neighbors ids
		/*!
		 * @param mesh The mesh
		 * @param mesh The vector with the vertexes ids
		 * @param mesh The set of the neighbor cells ids
		 */
		Cell2D( Geometry::Mesh2D * const mesh,
				const std::vector<UInt> & idVertexes,
				const std::set<UInt> & idNeighbors );
	
		//@}
	
		//! @name Get Methods
		//@{
		
		//! Get mesh (const)
		/*!
		 * @return The pointer to the related mesh
		 */
		Geometry::Mesh2D * const getMesh() const
			{ return M_mesh; }
	
		//! Get the vector containing the vertexes ids (const)
		/*!
		 * @return A constant reference to the vector of the vertexes ids
		 */
		const std::vector<UInt> & getVertexesVector() const
			{ return M_idVertexes; }
		
		//! Get the vector containing the vertexes ids
		/*!
		 * @return A reference to the vector of the vertexes ids
		 */
		std::vector<UInt> & getVertexesVector()
			{ return M_idVertexes; }
		
		//! Get the set containing the ids of the neighbor cells (const)
		/*!
		 * @return A constant reference to the set of the ids of the neighbor cells
		 */
		const std::set<UInt> & getNeighborsSet() const
			{ return M_idNeighbors; }
			
		//! Get the set containing the ids of the neighbor cells
		/*!
		 * @return A reference to the set of the ids of the neighbor cells
		 */
		std::set<UInt> & getNeighborsSet()
			{ return M_idNeighbors; }
	
		//! Get the cell centroid
		/*!
		 * @return A reference to the cell centroid
		 */
		const Geometry::Point2D & getCentroid() const
			{ return M_centroid; }
		
		//! Get the M_distanceCentroidFromFractureNetwork variable (const)
		/*!
		 * @return A constant reference to M_distanceCentroidFromFractureNetwork
		 */
		const Real & getDistanceFromFractureNetwork() const
			{ return M_distanceCentroidFromFractureNetwork; }
		
		//! Get the M_refinementOrder variable (const)
		/*!
		 * @return A constant reference to M_refinementOrder
		 */
		const UInt & getRefinementOrder() const
			{ return M_refinementOrder; }
			
		//@}
	
		//! @name Set Methods
		//@{
		
		//! Set the M_refinementOrder variable
		/*!
		 * @param refOrder The new value for M_refinementOrder
		 */	
		void setRefinementOrder( const UInt & refOrder )
			{ M_refinementOrder = refOrder; }

		//@}
	
		//! @name Methods
		//@{
		
		//! Insert a new vertex
		/*!
		 * @param idVertex The id of the new vertex
		 */
		void push_back(const UInt & idVertex)
			{ M_idVertexes.push_back(idVertex); }

		//! The number of vertexes of the cell
		/*!
		 * @return The number of vertexes
		 */
		UInt vertexesNumber() const
			{ return M_idVertexes.size(); }
	
		//! The number of neighbors of the cell
		/*!
		 * @return The number of neighbors
		 */
		UInt neighborsNumber() const
			{ return M_idNeighbors.size(); }

		//! The normalized vector normal to the edge defined by the two vertexes
		/*!
		 * @param idv1 The first vertex id
		 * @param idv2 The second vertex id
		 * @return The normalized vector normal to the edge (idv1,idv2)
		 */
		Geometry::Vector2D outerNormalToEdge( const UInt & idv1, const UInt & idv2 ) const;

		//! Test if a cell has a neighbor cell through the edge defined by the two vertexes
		/*!
		 * @param idv1 The first vertex id
		 * @param idv2 The second vertex id
		 * @param idv2 The variable in which the id of the neighbor cell is stored
		 * @return TRUE -> the neighbor cell exist
		 *		   FALSE -> the neighbor cell does not exist
		 */
		bool hasNeighborsThroughEdge( const UInt & vertex1, const UInt & vertex2,
								   	  UInt & idNeighbor) const;

		//! Computes the cell centroid and stores it in the M_centroid variable
		void computeCentroid();
	
		//! Computes the distance between the centroid and the fracture network fn
		/*!
		 * @param fn The fracture network
		 */
	/*	void computeCentroidDistanceFromFractureNetwork(const Geometry::FractureNetwork2D & fn)
			{ M_distanceCentroidFromFractureNetwork =
				  CGAL::sqrt( fn.squared_distance(M_centroid) ); }*/ //TODO tolto da me
		
		//! Computes the length of the longest edge of the cell
		/*!
		 * @param maxEdge A pair of two UInt storing the internal ids of the vertexes of the longest edge
		 * @return The length of the longest edge
		 */
		Real maxEdge( std::pair<UInt,UInt> & maxEdge ) const;
	
		//! Computes the signed area of the cell
		/*!
		 * @return The signed area of the cell
		 */
		Real signedArea() const;

		//! Computes the area of the cell
		/*!
		 * @return The area of the cell
		 */
		Real area() const;
		
		//! Computes the length of the diameter of the incircle
		/*!
		 * Precondition: cell must be a triangle
		 * @return The length of the diameter of the incircle
		 */
		Real incircleDiameter() const;

		//! Computes the length of the diameter of the circumcircle
		/*!
		 * Precondition: cell must be a triangle
		 * @return The length of the diameter of the circumcircle
		 */		
		Real circumcircleDiameter() const;
		
		//! Computes a quality index evaluating the cell shape
		/*!
		 * Precondition: cell must be a triangle
		 * @return The ratio circumcircleDiameter/incircleDiameter
		 */
		Real qualityIndex() const
			{ return circumcircleDiameter()/incircleDiameter(); }

		//! Test if a cell has a small angle
		/*!
		 * @param minAngle The cutoff value for the amplitude of the angles (in dregees)
		 * @return TRUE -> the cell has a small angle
		 * 		   FALSE -> the cell has a small angle
		 */
		bool has_small_angle( const Real & minAngle ) const;

		//! Test if a cell has a small angle and give a reference to the longest non fracture edge
		/*!
		 * @param minAngle The cutoff value for the amplitude of the angles (in dregees)
		 * @param edgeToBeFlipped A reference to the edge of the mesh that could be flipped
		 * @return TRUE -> the cell has a small angle
		 * 		   FALSE -> the cell has a small angle
		 */
		bool has_small_angle( const Real & minAngle,
							  Edge2D & edgeToBeFlipped ) const;
		
		//! Test if a cell is collapsed on one edge (the three vertexes are collinear)
		/*!
		 * @return TRUE -> the cell is degenerated
		 * 		   FALSE -> the cell is not degenerated
		 */
		bool is_degenerated_triangle() const;
		
		//! Display general information about the content of the class
		/*!
		 * @param out Specify the output format (std::cout by default)
		 */	
		void showMe(std::ostream & out=std::cout) const;
	
		//@}
		
	private:
	
		//! The pointer to the mesh containing this edge
		Geometry::Mesh2D * M_mesh;
		//! The vector containing the vertexes ids
		std::vector<UInt> M_idVertexes;
			// the vertexes have to be inserted in the right order: clockwise or counterclockwise
		//! The set storing the ids of the neighbor cells
		std::set<UInt> M_idNeighbors;
		//! The centroid
		Geometry::Point2D M_centroid;
		//! The distance between the centroid and a given fracture network
		Real M_distanceCentroidFromFractureNetwork;
		//! The refinement order (used only in cells of a CMesh)
		UInt M_refinementOrder;
	};
	

protected:
		/*!
			@class NodesConnectivity
    		This class implements the concept of nodes connectivity in a Mesh2D.
    	*/
	class NodesConnectivity {
	public:
		//! Define a connection of a mesh node
		typedef std::pair<UInt,Geometry::Mesh2D::Edge2D*> Connection;
		//! Define the map storing all the connections of a node
		typedef std::map<UInt,Geometry::Mesh2D::Edge2D*> ConnectionsMap;
		//! Define the set storing the cells sharing a node
		typedef std::set<UInt> NodePatch;
	
		//! Empty constructor
		NodesConnectivity();
		
		//! Constructor, getting the numeber of nodes in a mesh
		NodesConnectivity(const UInt & nGridNodes);
		
		//! Destructor
		~NodesConnectivity();
		
		//! Set the dimension of the vector (it is equal to the number of nodes in the mesh)
		/*!
		 * @param vectorSize The new size for the vector
		 */
		void setDimension(const UInt & vectorSize)
			{ M_vector.resize(vectorSize); }
		
		//! Clear the content of the vector (deallocate the memory)
		void clear();
		
		//! Get the map of the connection related to the given point (const)
		/*!
		 * @param idPoint The point id
		 * @return A constant reference to the map of the connection related to the given point
		 */
		const ConnectionsMap & getMapOfNeighbor( const UInt & idPoint ) const
			{ return *M_vector[idPoint]; }
			
		//! Get the map of the connection related to the given point
		/*!
		 * @param idPoint The point id
		 * @return A reference to the map of the connection related to the given point
		 */
		ConnectionsMap & getMapOfNeighbor( const UInt & idPoint )
			{ return *M_vector[idPoint]; }
		
		//! Get the set containing the ids of the cells sharing a point (const)
		/*!
		 * @param idPoint The point id
		 * @return A reference to the set containing the ids of the cells sharing that point
		 */
		const NodePatch & getNodePatch( const UInt & idPoint ) const
			{ return M_nodePatchsMap.at(idPoint); }
		
		//! Get the set containing the ids of the cells sharing a point
		/*!
		 * @param idPoint The point id
		 * @return A reference to the set containing the ids of the cells sharing that point
		 */
		NodePatch & getNodePatch( const UInt & idPoint )
			{ return M_nodePatchsMap.at(idPoint); }
		
		//! Insert a given connection with the connections of the given point
		/*!
		 * @param idPoint The point id
		 * @param conn The connection to be inserted
		 */
		void insert( const UInt & idPoint, const Connection & conn );
		
		//! Update the M_nodePatchsMap using the informations contained in M_vector
		void updateNodePatchsMap();
		
	private:
		//! The vector storing the ConnectionsMap of the mesh nodes
		std::vector<ConnectionsMap * > M_vector;
		std::map<UInt,NodePatch> M_nodePatchsMap;
	};

public:

	typedef std::vector<UInt> ApproxFracture;

	//! @name Constructor & Destructor
	//@{
	
	//! Empty constructor
	Mesh2D();

	//! Copy constructor
	Mesh2D(const Geometry::Mesh2D & mesh);
	
	//! Constructor, for a Cartesian mesh
	/*!
	 * @param Lx The length of the domain along the x direction
	 * @param Ly The length of the domain along the y direction
	 * @param Nx The number of elements along the x direction
	 * @param Ny The number of elements along the y direction
	 */
	Mesh2D(const Real & Lx, const Real & Ly, const UInt & Nx, const UInt & Ny);
		// Nx = # elementi direzione x (nodi in x = Nx+1)
		// Ny = # elementi direzione y (nodi in y = Ny+1)
	
	//! Constructor from Triangulation2D object
	Mesh2D(const Geometry::Triangulation2D & t);

	//! Destructor
	virtual ~Mesh2D();

	//@}
	
	//! @name Get Methods
	//@{
	
	//! Get the vector storing the nodes of the mesh (const)
	/*!
	 * @return A constant reference to the vector storing the nodes of the mesh
	 */
	const std::vector<Geometry::Point2D> & getNodesVector() const
		{ return M_nodes; }
//	std::vector<Geometry::Point2D> & getNodesVector()
//		{ return M_nodes; }
		
	//! Get the set storing the edges of the mesh (const)
	/*!
	 * @return A constant reference to the set storing the edges of the mesh
	 */
	const std::set<Edge2D> & getEdgesSet() const
		{ return M_edges; }
//	std::set<Geometry::Edge2D> & getEdgesSet()
//		{ return M_edges; }
		
	//! Get the map storing the cells of the mesh (const)
	/*!
	 * @return A constant reference to the map storing the cells of the mesh
	 */
	const std::map<UInt,Cell2D> & getCellsMap() const
		{ return M_cells; }

	//! Get the map storing the cells of the mesh
	/*!
	 * @return A reference to the map storing the cells of the mesh
	 */
	std::map<UInt,Cell2D> & getCellsMap()
		{ return M_cells; }
	
	//! Get thenodes connectivity of the mesh (const)
	/*!
	 * @return A constant reference to the nodes connectivity of the mesh
	 */
	const NodesConnectivity & getNodesConnectivity() const
		{ return M_nodesConnectivity; }
	
	const Geometry::FractureNetwork2D & getFn() const
		{ return M_fn; }
	
	const std::vector<ApproxFracture> & getApproxFn() const
		{ return M_approxFractureNetwork; }
	
	//@}
	
	//! @name Methods
	//@{
	
	//! Get the number of cells
	UInt Nelements() const
		{ return M_cells.size(); }
	
	//! Clear all the content of the class
	void clear()
		{ M_nodes.clear(); M_edges.clear(); M_cells.clear(); M_nodesConnectivity.clear(); }

	//! Compute the distance between each cell and a given fracture network
	/*!
	 * @param fn The fracture network
	 */
	void computeDistanceFieldFromFractureNetwork(const Geometry::FractureNetwork2D & fn);
		// VERIFY IF it could become a bottleneck for large problem!
	
	void addFractureNetwork( Geometry::FractureNetwork2D & fn )
		{ M_fn = fn; }
	
	//! Add to the mesh a new cell where all its vertexes are new nodes for the mesh
	/*!
	 * @param cellVertexes The vector containing the vertexes of the new cell
	 */
	void addDisjointCell(const std::vector<Geometry::Point2D> & cellVertexes);

	//! Add to the mesh a new cell where all its vertexes are new nodes for the mesh
	/*!
	 * @param cell The cell to be added
	 */
	void addDisjointCell(const Cell2D & cell);
	
	//! Add to the mesh a new cell where some of its vertexes can already be nodes of the mesh
	/*!
	 * @param cellVertexes The vector containing the vertexes of the new cell
	 */
	void addCell(const std::vector<Geometry::Point2D> & cellVertexes);
	
	//! Add to the mesh a new cell where some of its vertexes can already be nodes of the mesh
	/*!
	 * @param cell The cell to be added
	 */
	void addCell(const Cell2D & cell);
	
	//! Add to the mesh all the cell of the triangulation t
	/*!
	 * @param t The triangulation to be added
	 */
	void addTriangulation(const Geometry::Triangulation2D & t);
	
	//! Split the mesh into two sub-grids using a cutoff value on the distance field
	/*!
	 * @param inner The mesh containing all the cells with M_distanceCentroidFromFractureNetwork less than cutoff
	 * @param outer The mesh containing all the cells with M_distanceCentroidFromFractureNetwork greater than cutoff
	 * @param cutoff The cutoff parameter
	 */
	void splitMesh( Mesh2D & inner, Mesh2D & outer,
					const Real & cutoff ) const;
	
	//! Locate the position of a point into the mesh
	/*!
	 * @param point The point to be located
	 * @param idCell The variable where the method saves the id of the cell containing the point
	 * @return TRUE -> the point is into the mesh domain
	 * 		   FALSE -> the point is not into the mesh domain or an error occur with some cell
	 */
	bool findCellContainingPoint( const Geometry::Point2D & point, UInt & idCell ) const;
	
	void genericUnifyTriangle();
	
	//! Update the nodes connectivity object
	void buildNodesConnectivity();
	
	Real averageCellsQualityIndex() const;
	
	std::pair<UInt,Real> maximumCellsQualityIndex() const;
	
	bool checkCellsConsistency( std::vector<UInt> & degeneratedCellsIds ) const;
	
	//! Export the mesh in vtk format
	/*!
	 * (Use paraview to open the vtk file)
	 * @param filename The name of the vtk file created by this method
	 * @return TRUE -> operation ended correctly
				FALSE -> an error occurred
	 */
	bool exportVtk(const std::string & fileName) const;
	
	//! Export in vtk format the mesh with the distance field
	/*!
	 * (Use paraview to open the vtk file)
	 * @param filename The name of the vtk file created by this method
	 * @return TRUE -> operation ended correctly
				FALSE -> an error occurred
	 */
	bool exportCentroidsDistanceFromFractureNetworkVtk(const std::string & fileName) const;

	//! Export the fracture edges in vtk format
	/*!
	 * (Use paraview to open the vtk file)
	 * @param filename The name of the vtk file created by this method
	 * @return TRUE -> operation ended correctly
				FALSE -> an error occurred
	 */
	bool exportFractureEdgesVtk(const std::string & fileName) const;
	
	//! Export some cells of the mesh in vtk format
	/*!
	 * (Use paraview to open the vtk file)
	 * @param filename The name of the vtk file created by this method
	 * @return TRUE -> operation ended correctly
				FALSE -> an error occurred
	 */
	bool exportCellsVtk(const std::string & fileName, const std::vector<UInt> & idCells) const;
	
	bool exportApproxFractureNetworkUniqueVtk(const std::string & fileName) const;
	
	bool exportApproxFractureNetworkVtk(const std::string & fileName) const;
	
	bool exportApproxFractureVtk(const std::string & fileName, const UInt & f) const;
	
protected:
	//! The fracture network
	Geometry::FractureNetwork2D M_fn;
	//! The vector storing all the nodes of the mesh (with no repetition)
	std::vector<Geometry::Point2D> M_nodes;
	//! The set storing all the edges of the mesh (with no repetition)
	std::set<Edge2D> M_edges;
	//! The map storing all the cells of the mesh (with no repetition) 
	std::map<UInt,Cell2D> M_cells;
	//! The approx fracture network (for each fracture it contains the discretization nodes)
	std::vector<ApproxFracture> M_approxFractureNetwork;
	
	//! The object representing the nodes connectivity (updated statically be the user!)
	NodesConnectivity M_nodesConnectivity;
	
};

//! Overloading of operator< between two Edge2D objects (needed by the set M_edges)
/*!
 * @param e1 The first edge to be tested
 * @param e2 The second edge to be tested
 * @return TRUE -> e1 < e2
		   FALSE -> e1 >= e2
 */
bool operator<(const Geometry::Mesh2D::Edge2D & e1, const Geometry::Mesh2D::Edge2D & e2);

} // namespace Geometry

#endif /* MESH2D_HPP_ */
