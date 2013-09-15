 /*!
 *	@file Rigid_Mesh.hpp
 *	@brief Class for unstructured mesh.
 *
 *	@author Francesco Della Porta 
 *
 */ 
#ifndef RIGID_MESH_HPP_
#define RIGID_MESH_HPP_
 
#include "../src_Turconi/TypeDefinition.hpp"
#include "../src_Turconi/Mesh2D.hpp"
//#include "../src_Turconi/Mesh3D.hpp"
#include "Dimension.hpp"
#include "PolygonalDomain.hpp"
#include "../DP_Darcy/DarcyTypeDefinitions.hpp"

#include <cmath>
#include <set>
#include <vector>
#include <tuple>
#include <algorithm>
#include <map>

namespace Geometry{

/*!
	@class Rigid_Mesh
   	This class implements a container for Meshes of type Geometry::Mesh2D and Geometry::Mesh3D. Although those containers allow to make changes to the mesh easilly, they are not efficient in matrix building for solving PDE's. The Rigid_Mesh container doesn't permit any changes to the mesh but is very fast to discretize differential problems. Rigid_Mesh has a template parameter in order to decide at compilation time wheather to use methods for the 2D or the 3D case. This is done by choosing as template class T either Geometry::Dimension<2> or Geometry::DImension<3>. As Geometry::Mesh2D and Geometry::Mesh3D the Rigid_Mesh class is a container for mesh with arbitrary-shaped cells. It is furthermore, a container for domains with fractures.
*/
template <class T> 
class Rigid_Mesh{
public:

	/*!
		@typedef Generic_Facet
   		This typedefinition permits to treat Mesh2D::Edge2D and Mesh3D::Facet3D as a Generic_Facet.
   	*/
	typedef typename T::Generic_Facet Generic_Facet;
	/*!
		@typedef Generic_Point
   		This typedefinition permits to treat Geometry::Point2D and Geometry::Point3D as a Generic_Point.
   	*/
	typedef typename T::Generic_Point Generic_Point;
	/*!
		@typedef Generic_Vector
   		This typedefinition permits to treat Geometry::Vector2D and Geometry::Vector3D as a Generic_Facet.
   	*/
	typedef typename T::Generic_Vector Generic_Vector;
	/*!
		@typedef Generic_Facet
   		This typedefinition permits to treat Mesh2D and Mesh3D as a Generic_Mesh.
   	*/
	typedef typename T::Generic_Mesh Generic_Mesh;
	/*!
		@typedef Generic_Facet
   		This typedefinition permits to treat Mesh2D::Celll2D and Mesh3D::Cell3D as a Generic_Cell.
   	*/
	typedef typename T::Generic_Cell Generic_Cell;
	/*!
		@typedef Generic_Facet
   		This typedefinition permits to treat Mesh2D::Celll2D and Mesh3D::Cell3D as a Generic_Cell.
   	*/
//	typedef typename T::Generic_Point Point;
	/*!
		@typedef Generic_Facet
   		This typedefinition permits to treat 2D and 3D junctures of Fractures as Fracture_Juncture.
   	*/
	typedef typename T::Fracture_Juncture Fracture_Juncture;
	/*!
		@typedef Generic_Facet
   		This typedefinition permits to treat the container of the Facets of a Mesh2D and of a Mesh3D as Generic_FacetsContainer.
   	*/
	typedef std::set<Generic_Facet> Generic_FacetsContainer;
	/*!
		@typedef Point
   		This typedefinition permits to treat Generic_Point as a Point.
   	*/
	typedef Generic_Point Point;



	/*!
		@class Edge2D
   		This class implements the concept of Facet of a Rigid_Mesh.
		This class is dimension-dipendent just from the template parameter T.
   	*/
	class Facet{

	public:

		//! @name Constructor & Destructor
		//@{

		//! Constructor for a Facet from a facet of a Mesh2D or a Mesh3D.
		/*!
			@param generic_facet is a facet of a Mesh2D or a Mesh3D.
			@param mesh is a pointer to the mesh to which the Facet belongs
			@param old_to_new_map is a map from the old Id of the Neighbour-Cells in the Mesh2D or Mesh3D to the new Id of those cells in the Rigid_Mesh
			@param m_id is the Id of the facet in the Rigid_Mesh
		*/
		Facet (const Generic_Facet & generic_facet, Geometry::Rigid_Mesh<T> * const mesh, const std::map<UInt, UInt> &old_to_new_map, const UInt m_id);

		//! Copy constructor for a Facet from a Facet of another Rigid_Mesh
		Facet (const Facet& facet, Geometry::Rigid_Mesh<T> *const mesh);

		//! Copy constructor for a Facet from a Facet
		Facet (const Facet& facet);

		//! The empty constructor has been deleted
		Facet () = delete;

		//! Default Destructor
		~Facet () = default;

		//@}



		//! @name Get Methods
		//@{
		
		//! Get mesh (const)
		/*!
		 * @return A const pointer to the related mesh
		 */
		Geometry::Rigid_Mesh<T> *const getMesh () const
			{ return M_mesh; }

		//! Get vertexes Ids (const)
		/*!
		 * @return A reference to a vector of the Facet-vertexes-Ids
		 */
		const std::vector<UInt> & getVertexesIds () const
			{ return M_Vertexes_Ids; }

		//! Get center (const)
		/*!
		 * @return The baricenter of the Facet
		 */
		Generic_Point getCenter () const
			{ return M_center; }

		//! Get unsigned Normal (const)
		/*!
		 * @return The unsigned normal vector to the Facet
		 */
		Generic_Vector getUnsignedNormal () const
			{ return M_UnsignedNormal; }

		//! Get separated Cells Ids (const)
		/*!
		 * @return A vector reference with the Ids of the cells separated from the Facet
		 */
		const std::vector<UInt> & getSeparatedCellsIds () const
			{ return M_separatedCells_Ids; }

		//! Get size (const)
		/*!
		 * @return The size of the Facet: a lenght in 2D, an area in 3D
		 */
		Real size () const
			{ return M_size; }

		//! Get Id (const)
		/*!
		 * @return The Id of the Facet in the containing Rigid_Mesh
		 */
		UInt getId () const
			{ return M_Id; }
		//@}



		//! @name Methods
		//@{
		//! Display general information about the Facet
		/*!
		 * @param out Specify the output format (std::cout by default)
		 */	
		void showMe (std::ostream & out=std::cout) const;

		//@}

	protected:
		//! The pointer to the Rigid_Mesh containing the Facet
		Geometry::Rigid_Mesh<T> * M_mesh;
		//! The Facet Id in the containing Rigid_Mesh
		UInt M_Id;
		//! The vector of the Ids of the Facet's Vertexes
		std::vector<UInt> M_Vertexes_Ids;
		//! The vector of the Ids of the Cells separated by Facet
		std::vector<UInt> M_separatedCells_Ids;
		//! The size of the Facet: a lenght in 2D, an area in 3D
		Real M_size;
		//! The baricenter of the Facet
		Generic_Point M_center; 
		//! The unsigned normal-vector to the Facet
		Generic_Vector M_UnsignedNormal;

	protected:
		//! @name Protected-Methods
		//@{
		//! Computes the unsigned normal to the Facet in 3D
		/*!
		 * @param Is a type trait
		 */	
		void computeUnsignedNormal (const Dimension<3> );
		//! Computes the unsigned normal to the Facet in 2D
		/*!
		 * @param Is a type trait
		 */	
		void computeUnsignedNormal (const Dimension<2> );
		//! Is a function called by constructor in the 3D case
		/*!
		 * @param generic_facet is a facet of Mesh3D
		 * @param Is a type trait
		 */	
		void M_constructor (const Generic_Facet & generic_facet,
		 const Dimension<3> ); 
		//! Is a function called by constructor in the 2D case
		/*!
		 * @param generic_facet is a facet of Mesh2D or a Mesh3D
		 * @param Is a type trait
		 */	
		void M_constructor (const Generic_Facet & generic_facet,
		 const Dimension<2> ); 
		//! Is a function called by constructor which computes the Center of the Facet
		void computeCenter ();
		//@}
	};


	/*!
		@class Edge2D
   		This class implements the concept of Cell of a Rigid_Mesh.
		This class is dimension-dipendent just from the template parameter T.
   	*/
	class Cell{

	public:
		//! @name Constructor & Destructor
		//@{

		//! Constructor for a Cell from a cell of a Mesh2D or a Mesh3D.
		/*!
			@param generic_cell is a facet of a Mesh2D or a Mesh3D.
			@param mesh is a pointer to the mesh to which the Facet belongs
			@param m_id is the Id of the cell in the Rigid_Mesh
		*/
		Cell (const Generic_Cell & generic_cell, Geometry::Rigid_Mesh<T> *const mesh, const UInt m_id);

		//! Copy constructor for a Cell from a Cell of another Rigid_Mesh
		Cell (const Cell& cell,  Geometry::Rigid_Mesh<T> * const mesh);

		//! Copy constructor for a Cell from a Cell
		Cell (const Cell& cell);

		//! The empty constructor has been deleted
		Cell () = delete;

		//! Default Destructor
		~Cell() = default;

		//@}


		//! @name Get Methods
		//@{
		
		//! Get mesh (const)
		/*!
		 * @return A const pointer to the related mesh
		 */
		Geometry::Rigid_Mesh<T> *const getMesh () const
			{ return M_mesh; }

		//! Get Vertexes Ids (const)
		/*!
		 * @return A reference to a vector which contains the vertexes of the Cell
		 */
		const std::vector<UInt> & getVertexesIds () const
			{ return M_Vertexes_Ids; }

		//! Get Neighbors Ids (const)
		/*!
		 * @return A reference to a vector which contains the Ids of the neighour Cells
		 */
		const std::vector<UInt> & getNeighborsIds () const
			{ return M_Neighbors_Ids; }

		//! Get Centroid (const)
		/*!
		 * @return A reference to the centroid of the Cell
		 */
		const Generic_Point & getCentroid () const
			{ return M_centroid; }

		//! Get Vertexes number (const)
		/*!
		 * @return The number of the vertexes of the Cell
		 */
		UInt vertexesNumber () const
			{ return M_Vertexes_Ids.size (); }

		//! Get neighbours number (const)
		/*!
		 * @return The number of neighbours of the Cell
		 */
		UInt neighborsNumber() const
			{ return M_Neighbors_Ids.size (); }

		//! Get volume (const)
		/*!
		 * @return The volume of the Cell: a real volume in 3D, an area in 2D
		 */
		Real getVolume () const
			{ return M_volume; }

		//! Get cell Id (const)
		/*!
		 * @return The Id of the Cell in the containing Rigid_Mesh
		 */
		UInt getId () const
			{ return M_Id; }
		//@}



		//! @name Methods
		//@{
		//! Display general information about the Cell
		/*!
		 * @param out Specify the output format (std::cout by default)
		 */	
		void showMe (std::ostream & out=std::cout) const;

		//! Has Cell, idNeighbour as Neighbour through Facet facet_Id?
		/*!
		 * @return True if the Cell has idNeighbour as neighbour through Facet facet_Id, False if contrary 
		 * @param facet_Id is the Facet Id
		 * @param idNeighbour is the Cell Id
		 */	
		bool hasNeighborsThroughFacet (const UInt & facet_Id, const UInt & idNeighbor) const;
		//@}

		friend class Rigid_Mesh<T>;

	protected:
		//! The pointer to the Rigid_Mesh containing the Cell
		Geometry::Rigid_Mesh<T> * M_mesh;
		//! The Cell Id in the containing Rigid_Mesh
		UInt M_Id;
		//! The vector of the Ids of the Cell's Vertexes
		std::vector<UInt> M_Vertexes_Ids;
		//! The vector of the Ids of the Cells neighbours
		std::vector<UInt> M_Neighbors_Ids;
		//! The centroid of the Cell
		Generic_Point M_centroid;
		//! The size of the Cell: an area in 2D, a volume in 3D
		Real M_volume;

	protected:
		//! @name Protected-Methods
		//@{
		//! Is a function called by constructor in the 2D case
		/*!
		 * @param generic_cell is a cell of Mesh2D
		 * @param Is a type trait
		 */	
		void M_constructor (const Generic_Cell & generic_cell, const Dimension<2>); 
		//! Is a function called by constructor in the 3D case
		/*!
		 * @param generic_cell is a cell of Mesh3D
		 * @param Is a type trait
		 */	
		void M_constructor (const Generic_Cell & generic_cell, const Dimension<3>);
		//@}

	};



	/*!
		@class Facet_Id
   		This class is a container of a Facet-Id.
		This class is a base class.
   	*/
	class Facet_ID{

	protected:
		//! The Id of the referred Facet
		UInt Facet_Id;
		//! The pointer to the Rigid_Mesh containing the Facet
		Rigid_Mesh<T>* M_mesh;

	public:
		//! @name Constructor & Destructor
		//@{

		//! Constructor for a Facet_Id.
		/*!
			@param facet_id is the Id of a Facet of a Rigid_Mesh.
			@param mesh is a pointer to the mesh to which the Facet belongs
		*/
		Facet_ID (const UInt facet_id,  Geometry::Rigid_Mesh<T> *const mesh);

		//! The empty constructor has been deleted
		Facet_ID () = delete;

		//! void Destructor declared virtual
		virtual ~Facet_ID () {};
		//@}

		//! @name Get Methods
		//@{
		
		//! Get mesh (const)
		/*!
		 * @return A const pointer to the related mesh
		 */
		Geometry::Rigid_Mesh<T> *const getMesh () const
			{ return M_mesh; }
		
		//! Get Facet Id (const)
		/*!
		 * @return The const Facet-Id contained in the class
		 */
		const UInt getFacetId () const
			{return Facet_Id;}
		
		//! Get Facet (const)
		/*!
		 * @return A reference to a Facet whose Id is contained in the class
		 */
		const Rigid_Mesh<T>::Facet& getFacet () const
			{return M_mesh->getFacetsVector()[Facet_Id];}
		
		//! Get unsigned normal (const)
		/*!
		 * @return The unsigned normal vector to the Facet whose Id is contained in the class
		 */
		const Generic_Vector getUNormal () const
			{return M_mesh->getFacetsVector()[Facet_Id].getUnsignedNormal();}
		
		//! Get center (const)
		/*!
		 * @return The baricenter of the Facet whose Id is contained in the class
		 */
		const Generic_Point getCenter () const
			{return M_mesh->getFacetsVector()[Facet_Id].getCenter();}
		
		//! Get size (const)
		/*!
		 * @return The size of the Facet whose Id is contained in the class
		 */
		const Real getSize () const
			{return M_mesh->getFacetsVector()[Facet_Id].size();}
		
		//! Get separated cells (const)
		/*!
		 * @return The Id of the separated cell from Facet whose Id is contained in the class
		 */
		const std::vector<UInt>& getSeparated () const
			{return M_mesh->getFacetsVector()[Facet_Id].getSeparatedCellsIds();}
		//@}
	};


public:
	
	/*!
		@class Regular_Facet
   		This class is derived from the class Facet-Id.
		This class contains the Id of a Facet which is not a border-Facet and not a fracture-Facet.
   	*/
	class Regular_Facet : public Facet_ID {
	public:
		//! @name Constructor & Destructor
		//@{

		//! Constructor for a Regular_Facet given a Facet Id.
		/*!
			@param facet_id is the Id of a Facet of a Rigid_Mesh.
			@param mesh is a pointer to the mesh to which the Facet belongs
		*/
		Regular_Facet (const UInt facet_Id, Geometry::Rigid_Mesh<T> *const mesh);
		//! Copy constructor for a Regular_Facet given a regular_facet of another Rigid_Mesh.
		/*!
			@param regular_facet is a Regular_Facet of a different Rigid_Mesh.
			@param mesh is a pointer to the mesh to which the Facet belongs
		*/
		Regular_Facet (const Regular_Facet& regular_facet, Geometry::Rigid_Mesh<T> *const mesh);
		//! Copy constructor for a Regular_Facet given a regular_facet.
		Regular_Facet (const Regular_Facet &);
		//! Default destructor
		~Regular_Facet () = default;
	};
		//@}


	/*!
		@class Border_Facet
   		This class is derived from the class Facet-Id.
		This class contains the Id of a Facet which belongs to the border of the Domain.
   	*/
	class Border_Facet: public Facet_ID {
	protected:
		//! Identifies the Border
		UInt Border_Id;
	public:
		//! @name Constructor & Destructor
		//@{

		//! Constructor for a Border_Facet given a Facet Id and the relative border Id.
		/*!
			@param facet_id is the Id of a Facet of a Rigid_Mesh.
			@param border_id is the Id of the Domain-border to which the Facet belongs.
			@param mesh is a pointer to the mesh to which the Facet belongs
		*/
		Border_Facet (const UInt facet_Id, const UInt border_Id, Geometry::Rigid_Mesh<T> *const mesh);
		//! Copy constructor for a Border_Facet given a border_facet belonging to another Rigid_Mesh.
		Border_Facet (const Border_Facet& border_facet, /*const*/ Geometry::Rigid_Mesh<T> *const mesh);
		//! Copy constructor for a Border_Facet given a border_facet.
		Border_Facet (const Border_Facet &);
		//! Default destructor
		~Border_Facet () = default;
		//@}

		//! @name Get Methods
		//@{
		//! Get Border Id (const)
		/*!
		 * @return The Id of the border of the contained Facet
		 */
		UInt getBorderId() const
			{return Border_Id;}
		//@}
	};


	/*!
		@class Fracture_Facet
   		This class is derived from the class Facet-Id.
		This class contains the Id of a Facet which belongs to a fracture.
   	*/
	class Fracture_Facet: public Facet_ID {
	public:
		//! @name Constructor & Destructor
		//@{

		//! Constructor for a Fracture_Facet given a Facet Id and the Id of the represented fractures.
		/*!
			@param facet_Ids is a triple of long unsogned int: 1. Id of a Facet of a Rigid_Mesh; 2. Id as Fracture_Facet 3. Number of Cells in Rigid_Mesh  
			@param fracture_Ids is a vector with the Ids of the fractures represented by the Facet.
			@param mesh is a pointer to the mesh to which the Facet belongs
			@param generic_mesh is a const pointer to the Generic_Mesh, with the general informations about the fracture
		*/
		Fracture_Facet (const std::tuple<UInt,UInt,UInt> facet_Ids, const std::set<UInt>& fracture_Ids, Geometry::Rigid_Mesh<T> *const mesh, Generic_Mesh *const generic_mesh) ;
		//! Copy constructor for a Fracture_Facet given a fracture_facet belonging to another Rigid_Mesh.
		Fracture_Facet (const Fracture_Facet& fracture_facet, Geometry::Rigid_Mesh<T> *const mesh);
		//! Copy constructor for a Fracture_Facet
		Fracture_Facet (const Fracture_Facet&);
		//! Default destructor
		~Fracture_Facet () = default;
		//@}

		//! @name Get Methods
		//@{
		//! Permeability (const)
		/*!
		 * @return The permeability of the represented fractures
		 */
		const Real& Permeability () const
			{return M_permeability;}

		//! Aperture (const)
		/*!
		 * @return The Aperture of the represented fractures
		 */
		const Real& Aperture () const
			{return M_aperture;}

		//! Id as Fracture_Facet (const)
		/*!
		 * @return The Id as Fracture_Facet
		 */
		const UInt& getId () const
			{return M_Id;}

		//! Cells Number (const)
		/*!
		 * @return The number of cells in the containing Rigid_Mesh
		 */
		const UInt& getCellsnumber () const
			{return Cells_number;}

		//! Id as Cell (const)
		/*!
		 * @return The Id as Cell = Cells_Number + Id as Fracture_Facet
		 */
		const UInt getIdasCell () const
			{return (M_Id+Cells_number);}

		//! get Fractures Ids (const)
		/*!
		 * @return The Ids of the represented fractures
		 */
		const std::vector<UInt>& getFractureIds () const
			{return Fracture_Ids;}

		//! get Fracture Neighbours (const)
		/*!
		 * @return A map which associate to a border of the Facet the Ids of the Fracture_Facet neighbours from this border
		 */
		const std::map<Fracture_Juncture, std::vector<UInt> >& getFractureNeighbors () const
			{return Fracture_Neighbors;}
		//@}

		friend class Rigid_Mesh<T>; 

	private:
		//! Total number of Cell in the Rigid_Mesh containing this object
		UInt Cells_number;
		//! Id as Fracture_Facet of this Fracture_Facet 
		UInt M_Id;
		//! Permeability of the represented fractures
		Real M_permeability;
		//! Aperture of the represented fractures
		Real M_aperture;
		//! Ids of the represented fractures
		std::vector<UInt> Fracture_Ids;
		//! map which associate to a border of the Facet the Ids of the Fracture_Facet neighbours from this border
		std::map<Fracture_Juncture, std::vector<UInt> > Fracture_Neighbors;

	};


public:
	//! @name Constructor & Destructor
	//@{

	//! Constructor for a Rigid_Mesh given a Generic_Mesh and a Domain of type Geometry::Domain<T>.
	/*!
		@param generic_mesh is a reference to a Generic_Mesh, from which the Rigid_Mesh is constructed
		@param domain is Geometry::Domain<T> type and is necessary in order to give a tag to the border-facets
	*/
	Rigid_Mesh (Generic_Mesh & generic_mesh, Domain<T>& domain);
	//! Copy constructor
	Rigid_Mesh (const Rigid_Mesh &mymesh);
	//! empty constructor has been deleted
	Rigid_Mesh () = delete;
	//! Default destructor
	~Rigid_Mesh () = default;
	//@}


	//! @name Get Methods
	//@{

	//! Get nodes vector (const)
	/*!
	 * @return A reference to the vector containing the nodes of the mesh
	 */
	const std::vector<Point> & getNodesVector () const
		{ return M_nodes; }

	//! Get Facets vector (const)
	/*!
	 * @return A reference to the vector containing the Facets of the mesh
	 */
	const std::vector<Facet> & getFacetsVector () const
		{ return M_facets; }

	//! Get Cells vector (const)
	/*!
	 * @return A reference to the vector containing the Cells of the mesh
	 */
	const std::vector<Cell> & getCellsVector () const
		{ return M_cells; }

	//! Get internal Facets Ids vector (const)
	/*!
	 * @return A reference to the vector of the Regular_Facet
	 */
	const std::vector<Regular_Facet> & getInternalFacetsIdsVector () const
		{ return M_internalFacets; }

	//! Get border Facets Ids vector (const)
	/*!
	 * @return A reference to the vector of the Border_Facet
	 */
	const std::vector<Border_Facet> & getBorderFacetsIdsVector () const
		{ return M_borderFacets; }

	//! Get fracture Facets Ids vector (const)
	/*!
	 * @return A reference to the vector of the Fracture_Facet
	 */
	const std::vector<Fracture_Facet> & getFractureFacetsIdsVector () const
		{ return M_fractureFacets; }

	//! Get Domain (const)
	/*!
	 * @return A reference to the Domain of the mesh
	 */
	const Domain<T> & getDomain () const
		{ return M_domain; }

	//! Get size of biggest Facet (const)
	/*!
	 * @return The size of the biggest Facet
	 */
	const Real getMaxFacetSize () const
	{ 
		auto comp=[](Facet f1, Facet f2){return f1.size()<f2.size();};
		return std::max_element(M_facets.begin(),M_facets.end(), comp)->size(); 
	}

	//! Get size of smallest Facet (const)
	/*!
	 * @return The size of the smallest Facet
	 */
	const Real getMinFacetSize () const
	{ 
		auto comp=[](Facet f1, Facet f2){return f1.size()<f2.size();};
		return std::min_element(M_facets.begin(),M_facets.end(), comp)->size(); 
	}

	//! Get average size of Facet (const)
	/*!
	 * @return The average size of the Facets
	 */
	const Real getAveFacetSize () const
	{ 	
		Real ave = 0;
		auto sum = [&ave](Facet f1){ave+=f1.size();};
		for(auto it: M_facets) sum(it);
		return ave/M_facets.size();
	}
	//@}


	//! @name Methods
	//@{

	//! Display general information about the Rigid_Mesh
	/*!
	 * @param out Specify the output format (std::cout by default)
	 */	
	void showMe ( std::ostream & out=std::cout ) const;

	//! Saves the Rigid_Mesh in vtk format
	/*!
	 * @param filename Specify the name of the file where to save the mesh 
	 * @return True if the operation was successfull, False if not
	 */	
	bool exportVtk(const std::string & fileName) const;

	//! Appends in a file a solution of a PDE solved on this mesh in vtk format
	/*!
	 * @param sol Is the vector with the solution
	 * @param filename Specify the name of the file where to append the solution 
	 * @param label Specify the name of the variable which is saved ("pressure" by default)
	 * @param solType If "CELL" saves each component of the vector as it was the solution in a Cell, if "POINT" saves each component of the vector as it was the solution in a Point ("CELL" by default)
	 * @return True if the operation was successfull, False if not
	 */	
	bool appendSolutionToVtk (Darcy::Vector& sol, const std::string & fileName, const std::string & label = "pressure", const std::string & solType = "CELL") const;

	//! Converts Ids of Nodes to Nodes
	/*!
	 * @param pointsId A reference to a vector of nodes-Ids 
	 * @return A vector with the corresponding Points
	*/
	const std::vector<Generic_Point> IdToPoints (const std::vector<UInt>& pointsId);

	//@}


protected:
	//! @name Protected Methods
	//@{

	//! Is called by the constructor in the 2D case
	void M_constructor ( Generic_Mesh& generic_mesh, const Dimension<2> );
	//! Is called by the constructor in the 3D case
	void M_constructor ( Generic_Mesh& generic_mesh, const Dimension<3> );
	//! Shows the dimension of the Mesh in the 2D case
	/*!
	 * @param out Specify the output format (std::cout by default)
	 */	
	void showDimension ( const Dimension<2>, std::ostream & out=std::cout ) const;
	//! Shows the dimension of the Mesh in the 3D case
	/*!
	 * @param out Specify the output format (std::cout by default)
	 */	
	void showDimension ( const Dimension<3>, std::ostream & out=std::cout ) const;
	//! Prints the components of a point in the 2D case
	/*!
	 * @param generic_point Is the point whose components we want to print
	 * @param out Specify the output format (std::cout by default)
	 */	
	void showPoint ( const Geometry::Point2D& generic_point, std::ostream & out=std::cout ) const;
	//! Prints the components of a point in the 3D case
	/*!
	 * @param generic_point Is the point whose components we want to print
	 * @param out Specify the output format (std::cout by default)
	 */	
	void showPoint ( const Geometry::Point3D& generic_point, std::ostream & out=std::cout ) const;
	//! Builds the vector of Cells and is called by constructor
	/*!
	 * @param generic_mesh Is the generic_mesh we want to transform in a Rigid_Mesh
	 * @param old_to_new_map Is a map which binds the Id of a Cell in the originary Generic_Mesh to the Id in the Rigid_Mesh 
	 */	
	void CellsVectorBuilder ( Generic_Mesh & generic_mesh, std::map<UInt, UInt>& old_to_new_map );
	//! Builds the vector of Facets and is called by constructor in the 2D case
	/*!
	 * @param generic_mesh Is the generic_mesh we want to transform in a Rigid_Mesh
	 * @param old_to_new_map Is a map which binds the Id of a Cell in the originary Generic_Mesh to the Id in the Rigid_Mesh 
	 */	
	void FacetsVectorsBuilder ( Generic_Mesh & generic_mesh, const std::map<UInt, UInt>& old_to_new_map, const Dimension<2> );
	//! Builds the vector of Facets and is called by constructor in the 2D case
	/*!
	 * @param generic_mesh Is the generic_mesh we want to transform in a Rigid_Mesh
	 * @param old_to_new_map Is a map which binds the Id of a Cell in the originary Generic_Mesh to the Id in the Rigid_Mesh 
	*/
	void FacetsVectorsBuilder ( Generic_Mesh & generic_mesh, const std::map<UInt, UInt>& old_to_new_map, const Dimension<3> );
	//! Builds the vector of Facets and is called by constructor
	/*!
	 * @param old_to_new_map Is a map which binds the Id of a Cell in the originary Generic_Mesh to the Id in the Rigid_Mesh 
	 * @param facet_container Is the container of the facets in the generic_mesh we want to transform in a Rigid_Mesh
	*/
	void M_facetsVectorsBuilder(const std::map<UInt, UInt> & old_to_new_map, const Generic_FacetsContainer& facet_container);
	//! It is called by constructor and converts the Cells neighbours Id from old one to new one  
	/*!
	 * @param old_to_new_map Is a map which binds the Id of a Cell in the originary Generic_Mesh to the Id in the Rigid_Mesh 
	*/
	void AdjustCellNeighbors ( const std::map<UInt, UInt>& old_to_new_map );
	//! It is called by exportVtk and saves Cells on the Vtk in 2D 
	/*!
	 * @param filestr Specify the stream of the file where to save the mesh 
	 * @return True if the operation was successfull, False if not
	*/
	bool exportCellVtk (std::ostream & filestr, const Geometry::Dimension<2>) const;
	//! It is called by exportVtk and saves Cells on the Vtk in 3D 
	/*!
	 * @param filestr Specify the stream of the file where to save the mesh 
	 * @return True if the operation was successfull, False if not
	*/
	bool exportCellVtk (std::ostream & filestr, const Geometry::Dimension<3>) const;

protected:
	
	//! Vector of Nodes
	std::vector<Point> M_nodes;
	//! Vector of Facets
	std::vector<Facet> M_facets;
	//! Vector of Cells
	std::vector<Cell> M_cells;

	//! Vector of Regular_Facet
	std::vector<Regular_Facet> M_internalFacets;
	//! Vector of Border_Facet
	std::vector<Border_Facet> M_borderFacets;
	//! Vector of Fracture_Facet
	std::vector<Fracture_Facet> M_fractureFacets;

	//! Domain of the mesh
	Domain<T> M_domain;
};


// --------------------   Class Rigid_Mesh   --------------------
 
// ==================================================
// Protected Method
// ==================================================

	template<class T>
	void Rigid_Mesh<T>::showDimension (const Dimension<2>, std::ostream & out) const
	{
		out << "2D";
	}

	template<class T>
	void Rigid_Mesh<T>::showDimension (const Dimension<3>,  std::ostream & out) const
	{
		out << "3D";
	}

	template<class T>
	void Rigid_Mesh<T>::showPoint ( const Geometry::Point2D& generic_point, std::ostream & out) const
	{
		 out << generic_point.x() << " " << generic_point.y() << " 0" << std::endl;
	}

	template<class T>
	void Rigid_Mesh<T>::showPoint (const Geometry::Point3D& generic_point, std::ostream & out) const
	{
		 out << generic_point.x() << " " << generic_point.y() << " " << generic_point.z() << std::endl;
	}

	template<class T>
	void Rigid_Mesh<T>::CellsVectorBuilder (Generic_Mesh & generic_mesh,
	 std::map<UInt, UInt> & old_to_new_map)
	{
		UInt position = 0;
		for (auto it = generic_mesh.getCellsMap().begin();
		 it != generic_mesh.getCellsMap().end(); ++it)
			it->second.computeCentroid();
		for (auto it = generic_mesh.getCellsMap().begin();
		 it != generic_mesh.getCellsMap().end(); ++it)
		{
			old_to_new_map[it->first] = position;
			Cell cell(it->second, this, position);
			M_cells.push_back (cell);
			++ position;
		}	
	}
	

	template<class T>
	void Rigid_Mesh<T>::AdjustCellNeighbors ( const std::map<UInt, UInt> & old_to_new_map )
	{
		for (auto it = M_cells.begin(); it != M_cells.end(); ++it)
			for (auto Neighbors_it = it->M_Neighbors_Ids.begin();
			 Neighbors_it != it->M_Neighbors_Ids.end(); ++ Neighbors_it)
			{
				*(Neighbors_it) = old_to_new_map.at(*(Neighbors_it));
				#ifdef __MY__DEBUG__
				UInt comparison=-1;
				assert (*(Neighbors_it) != comparison);
				#endif
			}
	}
	

	template<class T>
	void Rigid_Mesh<T>::FacetsVectorsBuilder ( Generic_Mesh & generic_mesh,
	 const std::map<UInt, UInt> & old_to_new_map, const Dimension<2> )
	{
		M_facetsVectorsBuilder (old_to_new_map, generic_mesh.getEdgesSet());
		
		std::map<UInt, std::vector<UInt> > nodes_fracture_map;	

		for (auto it = M_fractureFacets.begin(); it != M_fractureFacets.end(); ++it)
			for (auto vertex_it = it->getFacet().getVertexesIds().begin(); vertex_it != it->getFacet().getVertexesIds().end(); ++vertex_it)		
				nodes_fracture_map[*(vertex_it)].emplace_back(it->getId());	

		for (auto it = M_fractureFacets.begin(); it != M_fractureFacets.end(); ++it)
			for (auto vertex_it = it->getFacet().getVertexesIds().begin();
				 vertex_it != it->getFacet().getVertexesIds().end(); ++vertex_it)	
			{	
				if (nodes_fracture_map.at(*(vertex_it)).size() > 1)
				{	
					it->Fracture_Neighbors[*(vertex_it)]=nodes_fracture_map.at(*(vertex_it));
					auto find_it = std::find(it->Fracture_Neighbors[*(vertex_it)].begin(),
						 it->Fracture_Neighbors[*(vertex_it)].end(), it->getId());
					it->Fracture_Neighbors[*(vertex_it)].erase(find_it);		
				}
			}
	}


	template<class T>	//ACHTUNG: STO ASSUMENDO CHE I PUNTI DI UNA FACET SIANO COLLEGATI DA UN SEGMENTO COL PUNTO PRIMA E QUELLO DOPO
	void Rigid_Mesh<T>::FacetsVectorsBuilder ( Generic_Mesh & generic_mesh,
	 const std::map<UInt, UInt> & old_to_new_map, const Dimension<3> )
	{
		M_facetsVectorsBuilder (old_to_new_map, generic_mesh.getFacetsSet());		

		std::map<Fracture_Juncture, std::vector<UInt> > nodes_fracture_map;	

		for (auto it = M_fractureFacets.begin(); it != M_fractureFacets.end(); ++it)
		{
			auto vertex_it2 = it->getFacet().getVertexesIds().rbegin();
			for (auto vertex_it = it->getFacet().getVertexesIds().begin();
				 vertex_it != it->getFacet().getVertexesIds().end(); ++vertex_it)
			{		
				//The if state is necessary in order to order the vertexes!!
				if (*(vertex_it2) < *(vertex_it))
					nodes_fracture_map[(*(vertex_it2), *(vertex_it))].emplace_back(it->getId());
				else			
					nodes_fracture_map[(*(vertex_it), *(vertex_it2))].emplace_back(it->getId());
				vertex_it2 = vertex_it;	
			}	
		}

		for (auto it = M_fractureFacets.begin(); it != M_fractureFacets.end(); ++it)
		{
			auto vertex_it2 = it->getFacet().getVertexesIds().rbegin();
			for (auto vertex_it = it->getFacet().getVertexesIds().begin();
				 vertex_it != it->getFacet().getVertexesIds().end(); ++vertex_it)	
			{	
				Fracture_Juncture m_junct;
				if (vertex_it < vertex_it2)
					m_junct = std::make_pair(*(vertex_it), *(vertex_it2));
				else
					m_junct = std::make_pair(*(vertex_it2), *(vertex_it));
				if (nodes_fracture_map.at(m_junct).size() > 1)
				{
					it->Fracture_Neighbors[m_junct] = nodes_fracture_map.at(m_junct);
					auto find_it = std::find(it->Fracture_Neighbors[m_junct].begin(),
					 it->Fracture_Neighbors[m_junct].end(), it->getId());
					it->Fracture_Neighbors[m_junct].erase(find_it);
				}
				vertex_it2 = vertex_it;	
			}
		}
	}


	template<class T>
	void Rigid_Mesh<T>::M_facetsVectorsBuilder (const std::map<UInt,UInt> & old_to_new_map, const Generic_FacetsContainer& facet_container)
	{
		UInt position = 0;
		UInt fractureId = 0;
	
		for(auto it = facet_container.begin();
		 it != facet_container.end(); ++it)
		{
			Facet facet(*(it), this, old_to_new_map, position);
			M_facets.push_back(facet);

			if (it->isFracture())
			{
				std::tuple<UInt,UInt,UInt> fracture_Ids;
				std::get<0>(fracture_Ids) = position;
				std::get<1>(fracture_Ids) = fractureId;
				std::get<2>(fracture_Ids) = M_cells.size();

				Fracture_Facet fracture_facet ( fracture_Ids,
					 it->getRepresentedFractureIds(), this, it->getMesh()); 
				M_fractureFacets.push_back(fracture_facet);
				++fractureId;				
			}
			else
			{
				if (it->getSeparatedCells().size() == 1)
				{
					Border_Facet border_facet (position, M_domain.getBorderId(this->IdToPoints(M_facets.rbegin()->getVertexesIds())), this);
					UInt m_check = -1;
					assert(border_facet.getBorderId() != m_check);
					M_borderFacets.push_back (border_facet);
				}
				else
				{
					#ifdef __MY__DEBUG__
					assert (it->getSeparatedCells().size() == 2);
					#endif
	
					M_internalFacets.emplace_back(Regular_Facet (position, this));
				}	
			}
			++ position;		
		}			
		#ifdef __MY__DEBUG__
		assert (M_borderFacets.size()+M_fractureFacets.size()+M_internalFacets.size() == M_facets.size());
		#endif
	}	


	template<class T>
	bool Rigid_Mesh<T>::exportCellVtk (std::ostream & filestr, const Geometry::Dimension<2>) const
	{
		for ( auto cells_it : M_cells)
		{
			switch (cells_it.vertexesNumber())
			{
				case 3:
					filestr << "5" << std::endl;
					break;
				case 4:
					filestr << "9" << std::endl;
					break;
				case 5:
				case 6:
				case 7:
				case 8:
					filestr << "7" << std::endl;
					break;
				default:
					std::cerr << "   *** Error *** : Rigid_Mesh in 2D...exportVtk() failed, ";
					std::cerr << "unexpected cell geometry";
					return 0;
			}
		}
		for ( auto fractures_it : M_fractureFacets)
		{
			if (fractures_it.getFacet().getVertexesIds().size() == 2)
			{
				filestr << "3" << std::endl;
			}
			else
			{
					std::cerr << "   *** Error *** : Rigid_Mesh in 2D...exportVtk() failed, ";
				std::cerr << "unexpected cell geometry";
				return 0;
			}
		}
		return 1;
	}


	template<class T>
	bool Rigid_Mesh<T>::exportCellVtk (std::ostream & filestr , const Geometry::Dimension<3>) const
	{
		for ( auto cells_it : M_cells)
		{
			switch (cells_it.vertexesNumber()){
				case 4:
					filestr << "10" << std::endl;
					break;
				case 5:
					filestr << "14" << std::endl;
					break;
				case 6:
					filestr << "13" << std::endl;
					break;
				case 8:
					filestr << "12" << std::endl;
					break;
				default:
					std::cerr << "   *** Error *** : Rigid_Mesh in 3D...exportVtk() failed, ";
					std::cerr << "unexpected cell geometry";
				return 0;
			}
		}
		for ( auto fractures_it : M_fractureFacets)
		{
			switch (fractures_it.getFacet().getVertexesIds().size()){
				case 3:
					filestr << "5" << std::endl;
					break;
				case 4:
					filestr << "9" << std::endl;
					break;
				case 5:
				case 6:
				case 7:
				case 8:
					filestr << "7" << std::endl;
					break;
				default:
					std::cerr << "   *** Error *** : Rigid_Mesh in 3D...exportVtk() failed, ";
					std::cerr << "unexpected cell geometry";
				return 0;
			}
		}
		return 1;
	}


// ==================================================
// Constructors & Destructor
// ==================================================	

	template<class T>
	Rigid_Mesh<T>::Rigid_Mesh (Generic_Mesh & generic_mesh, Domain<T>& domain):
	 M_nodes (generic_mesh.getNodesVector()), M_domain(domain)
	{
		std::map<UInt, UInt> old_to_new_map;
		CellsVectorBuilder (generic_mesh, old_to_new_map);
		AdjustCellNeighbors (old_to_new_map);
		FacetsVectorsBuilder (generic_mesh, old_to_new_map, T());
	}

		
	template<class T>
	Rigid_Mesh<T>::Rigid_Mesh(const Rigid_Mesh<T> &mymesh):M_nodes(mymesh.getNodesVector()), M_domain(mymesh.getDomain())
	{
		for (auto it = mymesh.getCellsVector().begin(); it != mymesh.getCellsVector().end(); ++it)
		{
			Cell my_cell (*(it), this);
			M_cells.push_back (my_cell);
		}

		for (auto it = mymesh.getFacetsVector().begin(); it != mymesh.getFacetsVector().end(); ++it)
		{
			Facet my_facet (*(it), this);
			M_facets.push_back (my_facet);
		}

		for (auto it = mymesh.getBorderFacetsIdsVector().begin(); it != mymesh.getBorderFacetsIdsVector().end(); ++it)
		{
			M_borderFacets.emplace_back (Border_Facet (*(it), this));
		}

		for (auto it = mymesh.getInternalFacetsIdsVector().begin(); it != mymesh.getInternalFacetsIdsVector().end(); ++it)
		{
			M_internalFacets.emplace_back (Regular_Facet (*(it), this));
		}

		for (auto it = mymesh.getFractureFacetsIdsVector().begin(); it != mymesh.getFractureFacetsIdsVector().end(); ++it)
		{
			M_fractureFacets.emplace_back (Fracture_Facet (*(it), this));
		}
	}

	
// ==================================================
// Methods
// ==================================================	

	template<class T>
	bool Rigid_Mesh<T>::exportVtk(const std::string & fileName) const
	{
		std::fstream filestr;
		
		filestr.open (fileName.c_str(), std::ios_base::out);
	
		if (filestr.is_open())
		{
			std::cout << std::endl << " File: " << fileName << ", successfully opened";
		}
		else
		{
			std::cerr << std::endl
					  << " *** Error: file not opened *** " << std::endl << std::endl;
			return  0;
		}
		
		std::cout << std::endl << " Exporting Rigid Mesh in Vtk format... " << std::endl;
		std::cout << " The mesh is in " ;
		showDimension(T());
		std::cout << std::endl;
		
		UInt nCells = M_cells.size() + M_fractureFacets.size();
		UInt nPoints = M_nodes.size();
		
		// Header
		filestr << "# vtk DataFile Version 3.1" << std::endl;
		filestr << "this is a file created for Paraview" << std::endl;
		filestr << "ASCII" << std::endl;
		filestr << "DATASET UNSTRUCTURED_GRID" << std::endl;
		filestr << std::endl;
		
		filestr << std::scientific << std::setprecision(10);
		
		// Pointdata
		filestr << "POINTS " << nPoints << " double" << std::endl;
		for( auto point_it : M_nodes)
		{
			showPoint ( point_it, filestr );
		}
		filestr << std::endl;
			
		// Celldata
		UInt nVal = 0;
		for( auto cells_it : M_cells)
			nVal += cells_it.vertexesNumber() +1;
		for( auto fractures_it : M_fractureFacets)
			nVal += fractures_it.getFacet().getVertexesIds().size() +1;
		
		filestr << "CELLS " << nCells << " " << nVal << std::endl;
		
		for( auto cells_it : M_cells)
		{
			filestr << cells_it.vertexesNumber();
			for( auto nodes_it : cells_it.getVertexesIds())
			{
				filestr << " " << nodes_it;
			}
			filestr << std::endl;
		}
		for( auto fractures_it : M_fractureFacets)
		{
			filestr << fractures_it.getFacet().getVertexesIds().size();
			for( auto nodes_it : fractures_it.getFacet().getVertexesIds())
			{
				filestr << " " << nodes_it;
			}
			filestr << std::endl;
		}
		filestr << std::endl;
		
		filestr << "CELL_TYPES " << nCells << std::endl;
		bool flag = exportCellVtk(filestr, T());		
		filestr.close();
		
		return flag;
	}


	template<class T>
	bool Rigid_Mesh<T>::appendSolutionToVtk (Darcy::Vector& sol, const std::string & fileName, const std::string & label, const std::string & solType) const
	{
		std::fstream filestr;
		
		filestr.open (fileName.c_str(), std::ios_base::out | std::ios_base::app	);
	
		if (filestr.is_open())
		{
			std::cout << std::endl << " File: " << fileName << ", successfully opened";
		}
		else
		{
			std::cerr << std::endl
					  << " *** Error: file not opened *** " << std::endl << std::endl;
			return  0;
		}
		
		std::cout << std::endl << " Appending solution to Vtk file... " << std::endl << std::endl;
		
		UInt N = sol.size();
				
		// Solution
		filestr << std::endl;
		if (solType != "CELL" && solType != "POINT")
			std::cerr << std::endl  << " *** Error: wrong setting of solType *** " << std::endl << std::endl;
		filestr << solType <<"_DATA " << N << std::endl;
		filestr << "SCALARS " << label << " double "<< 1 << std::endl;
		filestr	<< "LOOKUP_TABLE" << " " << "default" << std::endl;	
		for( UInt sol_it = 0; sol_it < N; ++sol_it )
		{
			filestr << sol(sol_it) << std::endl;
		}
		filestr.close();
		return 1;
	}


	template<class T>
	void Rigid_Mesh<T>::showMe ( std::ostream & out ) const
	{
		out << "Type = Mesh";
		this->showDimension(T(), out);
		out <<" : "<< std::endl;
		out << "code = " << this << std::endl;
		out << "Number of Nodes: " << M_nodes.size() <<std::endl;
		out << "Number of Cells: " << M_cells.size() <<std::endl;
		out << "Number of Facets: " << M_facets.size() <<std::endl;
		out << "Number of Border-Facets: " << M_borderFacets.size() <<std::endl;
		out << "Number of Fracture-Facets: " << M_fractureFacets.size() <<std::endl;
		out << "Number of Standard-Facets: " << M_internalFacets.size() <<std::endl;

	}

	template<class T>
	const std::vector<typename T::Generic_Point> Rigid_Mesh<T>::IdToPoints (const std::vector<UInt>& pointsIds)
	{
		std::vector<Generic_Point> points;
		for (auto iter : pointsIds)
			points.emplace_back(M_nodes[iter]);
		return points;
	}













// --------------------   Class Facet   --------------------
 
// ==================================================
// Constructors & Destructor
// ==================================================

	template<class T>
	Rigid_Mesh<T>::Facet::Facet(const Generic_Facet & generic_facet, Geometry::Rigid_Mesh<T> *const mesh, const std::map<UInt,UInt> &old_to_new_map, const UInt m_id): M_mesh(mesh), M_Id(m_id)
	{
		M_constructor(generic_facet, T() );
		for(auto iterator = generic_facet.getSeparatedCells().begin(); iterator != generic_facet.getSeparatedCells().end(); 
		++iterator)
		{
			M_separatedCells_Ids.emplace_back(old_to_new_map.at(*(iterator)));
		}
		computeCenter();
		if (mesh == NULL)
			std::cerr << "NULL PTR GIVEN" << std::endl;
		computeUnsignedNormal(T());
	}


	template<class T>
	Rigid_Mesh<T>::Facet::Facet(const Facet& facet): M_mesh(facet.getMesh()), M_Id(facet.getId()), M_Vertexes_Ids(facet.getVertexesIds()), 
		M_separatedCells_Ids(facet.getSeparatedCellsIds()), M_size(facet.size()), M_center(facet.getCenter()),
		M_UnsignedNormal(facet.getUnsignedNormal()){}


	template<class T>
	Rigid_Mesh<T>::Facet::Facet(const Facet& facet, Geometry::Rigid_Mesh<T> *const mesh): M_mesh(mesh), M_Id(facet.getId()), M_Vertexes_Ids(facet.getVertexesIds()), 
		M_separatedCells_Ids(facet.getSeparatedCellsIds()), M_size(facet.size()), M_center(facet.getCenter()),
		M_UnsignedNormal(facet.getUnsignedNormal()){}


// ==================================================
// Protected Methods
// ==================================================

	template<class T>
	void Rigid_Mesh<T>::Facet::computeCenter()
	{
		int N = M_Vertexes_Ids.size();
		Generic_Point sum = CGAL::ORIGIN;
		for ( int j = 0; j < N; ++j)
			sum = sum + (M_mesh->getNodesVector()[M_Vertexes_Ids[j]]-CGAL::ORIGIN)/N;
		this->M_center =  sum;
	}


	template<class T>
	void Rigid_Mesh<T>::Facet::computeUnsignedNormal( const Dimension<2> )
	{
		Generic_Point a;
		Generic_Point b;
		a = M_mesh->getNodesVector()[M_Vertexes_Ids[0]];
		b = M_mesh->getNodesVector()[M_Vertexes_Ids[1]];
		Generic_Vector tangent;
		tangent = b - a;
		Generic_Vector normal (tangent.y(), -tangent.x());
		Real norm = sqrt(normal.x()*normal.x()+normal.y()*normal.y());//si noti che si usa il per e non l'elevamento
		normal=normal/norm;
		#ifdef __MY__DEBUG__
		assert(normal.x()*tangent.x() + normal.y()*tangent.y() < 1.e-14);
		assert(normal.x()*normal.x() + normal.y()*normal.y() < 1. + 1.e-14);
		assert(normal.x()*normal.x() + normal.y()*normal.y() > 1. - 1.e-14);
		#endif
		this->M_UnsignedNormal=normal;
	}


	template<class T>
	void Rigid_Mesh<T>::Facet::computeUnsignedNormal( const Dimension<3> )
	{
		Generic_Vector tangent_1, tangent_2;
		tangent_1 = M_mesh->getNodesVector()[M_Vertexes_Ids[1]] - M_mesh->getNodesVector()[M_Vertexes_Ids[0]];
		tangent_2 = M_mesh->getNodesVector()[M_Vertexes_Ids[2]] - M_mesh->getNodesVector()[M_Vertexes_Ids[0]];

		Generic_Vector normal( tangent_1.y()*tangent_2.z() - tangent_1.z()*tangent_2.y(),
								tangent_1.z()*tangent_2.x() - tangent_1.x()*tangent_2.z(),
								tangent_1.x()*tangent_2.y() - tangent_1.y()*tangent_2.x());

		Real norm = sqrt(normal.x()*normal.x() + normal.y()*normal.y() + normal.z()*normal.z());
		normal=normal/norm;

		#ifdef __MY__DEBUG__
		assert(normal.x()*tangent_1.x() + normal.y()*tangent_1.y() + normal.z()*tangent_1.z() < 1.e-16);
		assert(normal.x()*tangent_2.x() + normal.y()*tangent_2.y() + normal.z()*tangent_2.z() < 1.e-16);
		assert(normal.x()*normal.x() + normal.y()*normal.y() + normal.z()*normal.z() < 1. + 1.e-14);
		assert(normal.x()*normal.x() + normal.y()*normal.y() + normal.z()*normal.z() > 1. - 1.e-14);
		#endif
		this->M_UnsignedNormal=normal;
	}


	template<class T>
	void Rigid_Mesh<T>::Facet::M_constructor(const Generic_Facet & generic_facet, const Dimension<3> )
	{
		M_Vertexes_Ids=generic_facet.getVertexesVector();
		M_size = generic_facet.area();
	}


	template<class T>
	void Rigid_Mesh<T>::Facet::M_constructor(const Generic_Facet & generic_facet, const Dimension<2> ) 
	{
		M_Vertexes_Ids.emplace_back(generic_facet.getP1());
		M_Vertexes_Ids.emplace_back(generic_facet.getP2());
		M_size = generic_facet.length();
	}


// ==================================================
// Methods
// ==================================================

	template<class T>
	void Rigid_Mesh<T>::Facet::showMe(std::ostream  & out) const
	{
		out << "Type = Facet";
		M_mesh->showDimension(T(), out); 
		out << " :" << std::endl;
		out << "ID: " << M_Id << std::endl;
		out << "M_mesh = " << M_mesh << std::endl;
		out << "Number of Nodes: " << M_Vertexes_Ids.size()<<" : "<<std::endl;
		UInt counter=1;
		for(auto j = M_Vertexes_Ids.begin(); j!=M_Vertexes_Ids.end(); ++j)
		{
			out << "M_idP"<<counter<<" = " <<*(j)<<"\t";
			++counter;
		}
		out << std::endl;
		out << "Size: " << M_size << std::endl;
		out << "Center: ";
		M_mesh->showPoint(M_center, out);
		out << "M_separatedCells : size = " << M_separatedCells_Ids.size();
		out << "    [ ";
		if( !M_separatedCells_Ids.empty() )
			for( auto it = M_separatedCells_Ids.begin();
				 it != M_separatedCells_Ids.end(); ++it )
			out << *it << " ";
		out << "] " << std::endl;
	}









	
// --------------------   Class Cell   --------------------
 
// ==================================================
// Constructors & Destructor
// ==================================================

	template<class T>
	Rigid_Mesh<T>::Cell::Cell( const Generic_Cell & generic_cell, Geometry::Rigid_Mesh<T>*const mesh, const UInt m_id): M_mesh(mesh), M_Id(m_id)
		 , M_Vertexes_Ids(generic_cell.getVertexesVector())
	{
		for(auto it = generic_cell.getNeighborsSet().begin(); it != generic_cell.getNeighborsSet().end(); ++it)
			M_Neighbors_Ids.emplace_back(*it);
		M_centroid = generic_cell.getCentroid();
		M_constructor(generic_cell, T());
	}


	template<class T>
	Rigid_Mesh<T>::Cell::Cell(const Cell& cell, Geometry::Rigid_Mesh<T> *const mesh): M_mesh(mesh), M_Id(cell.getId()), M_Vertexes_Ids(cell.getVertexesIds()), 
		M_Neighbors_Ids(cell.getNeighborsIds()), M_centroid(cell.getCentroid()), M_volume(cell.getVolume())
	{}

	template<class T>
	Rigid_Mesh<T>::Cell::Cell(const Cell& cell): M_mesh(cell.getMesh()), M_Id(cell.getId()), M_Vertexes_Ids(cell.getVertexesIds()), 
		M_Neighbors_Ids(cell.getNeighborsIds()), M_centroid(cell.getCentroid()), M_volume(cell.getVolume())
	{}


// ==================================================
// Protected Methods
// ==================================================

	template<class T>
	void Rigid_Mesh<T>::Cell::M_constructor(const Generic_Cell & generic_cell, const Dimension<3> )
	{
		M_volume=generic_cell.volume();
	}


	template<class T>
	void Rigid_Mesh<T>::Cell::M_constructor(const Generic_Cell & generic_cell, const Dimension<2> ) 
	{
		M_volume=generic_cell.area();
	}


// ==================================================
// Methods
// ==================================================

	template<class T>
	bool Rigid_Mesh<T>::Cell::hasNeighborsThroughFacet( const UInt & facet_Id, const UInt & idNeighbor) const
	{
		if(M_mesh->getFacetsVector()[facet_Id].getSeparatedCellsIds()[0] == idNeighbor
			|| M_mesh->getFacetsVector()[facet_Id].getSeparatedCellsIds()[1] == idNeighbor)
		{
			return true;
		}
		return false;
	}


	template<class T>
	void Rigid_Mesh<T>::Cell::showMe(std::ostream  & out) const
	{
		out << "Type = Cell";
		M_mesh->showDimension(T(), out); 
		out << " :" << std::endl;;
		out << "ID: " << M_Id << std::endl;
		out << "M_mesh = " << M_mesh << std::endl;
		out << "Number of Nodes: " << M_Vertexes_Ids.size()<<" : "<<std::endl;
		UInt counter=1;
		for(auto j = M_Vertexes_Ids.begin(); j!=M_Vertexes_Ids.end(); ++j)
		{
			out << "M_idP"<<counter<<" = " <<*(j)<<"\t";
			++counter;
		}
		out << std::endl;
		out << "Volume: " << M_volume << std::endl;
		out << "Centroid: ";
		M_mesh->showPoint(M_centroid, out);
		out << "M_Neighbors_Ids : size = " << M_Neighbors_Ids.size();
		out << "    [ ";
		if( !M_Neighbors_Ids.empty() )
			for( auto it = M_Neighbors_Ids.begin();
				 it != M_Neighbors_Ids.end(); ++it )
			out << *it << " ";
		out << "] " << std::endl;
	}





// --------------------   Class Facet_ID   --------------------
 
// ==================================================
// Constructors & Destructor
// ==================================================

	template<class T>
	Rigid_Mesh<T>::Facet_ID::Facet_ID(const UInt facet_Id, Geometry::Rigid_Mesh<T> *const mesh):
	  Facet_Id(facet_Id), M_mesh(mesh){}


// --------------------   Class Regular_Facet   --------------------
 
// ==================================================
// Constructors & Destructor
// ==================================================

	template<class T>
	Rigid_Mesh<T>::Regular_Facet::Regular_Facet(const UInt facet_Id, Geometry::Rigid_Mesh<T> *const mesh):
	  Facet_ID(facet_Id, mesh){}


	template<class T>
	Rigid_Mesh<T>::Regular_Facet::Regular_Facet(const Regular_Facet& regular_facet):
	  Facet_ID(regular_facet.getFacetId(), regular_facet.getMesh())
	{}


	template<class T>
	Rigid_Mesh<T>::Regular_Facet::Regular_Facet(const Regular_Facet& regular_facet, Geometry::Rigid_Mesh<T> *const mesh):
	  Facet_ID(regular_facet.getFacetId(), mesh)
	{}



// --------------------   Class Border_Facet   --------------------
 
// ==================================================
// Constructors & Destructor
// ==================================================

	template <class T>
	Rigid_Mesh<T>::Border_Facet::Border_Facet(const UInt facet_Id, const UInt border_Id, Geometry::Rigid_Mesh<T> *const mesh):  
	  Facet_ID(facet_Id, mesh),	Border_Id(border_Id){}

	template <class T>
	Rigid_Mesh<T>::Border_Facet::Border_Facet(const Border_Facet& border_facet, Geometry::Rigid_Mesh<T> *const mesh):
	  Facet_ID(border_facet.getFacetId(), mesh), Border_Id(border_facet.getBorderId()){}

	template <class T>
	Rigid_Mesh<T>::Border_Facet::Border_Facet(const Border_Facet& border_facet):
	  Facet_ID(border_facet.getFacetId(), border_facet.getMesh()), Border_Id(border_facet.getBorderId()) {}


	
// --------------------   Class Fracture_Facet   --------------------
 
// ==================================================
// Constructors & Destructor
// ==================================================

	template <class T>
	Rigid_Mesh<T>::Fracture_Facet::Fracture_Facet(const Fracture_Facet& fracture_facet, Geometry::Rigid_Mesh<T> *const mesh):
	  Facet_ID(fracture_facet.getFacetId(), mesh), Cells_number(fracture_facet.getCellsnumber()), M_Id(fracture_facet.getId()), M_permeability(fracture_facet.Permeability()), M_aperture(fracture_facet.Aperture()), Fracture_Ids(fracture_facet.getFractureIds()), Fracture_Neighbors(fracture_facet.getFractureNeighbors()){}


	template <class T>
	Rigid_Mesh<T>::Fracture_Facet::Fracture_Facet(const Fracture_Facet& fracture_facet):
	  Facet_ID(fracture_facet.getFacetId(), fracture_facet.getMesh()), Cells_number(fracture_facet.getCellsnumber()), M_Id(fracture_facet.getId()), M_permeability(fracture_facet.Permeability()), M_aperture(fracture_facet.Aperture()), Fracture_Ids(fracture_facet.getFractureIds()), Fracture_Neighbors(fracture_facet.getFractureNeighbors()){}


	template <class T>
	Rigid_Mesh<T>::Fracture_Facet::Fracture_Facet (const std::tuple<UInt,UInt,UInt> facet_Ids, const std::set<UInt>& fracture_Ids, Geometry::Rigid_Mesh<T> *const mesh, Generic_Mesh *const generic_mesh):
	 Facet_ID(std::get<0>(facet_Ids), mesh), Cells_number(std::get<2>(facet_Ids)), M_Id(std::get<1>(facet_Ids))
	{
		for (auto it = fracture_Ids.begin(); it != fracture_Ids.end(); ++it)
			Fracture_Ids.emplace_back(*(it));

		M_aperture = 0;
		M_permeability = 0;

		for (auto it = Fracture_Ids.begin(); it != Fracture_Ids.end(); ++it)
		{
			M_aperture += generic_mesh->getFn().getNetwork()[*(it)].aperture();
			M_permeability +=  generic_mesh->getFn().getNetwork()[*(it)].aperture()*generic_mesh->getFn().getNetwork()[*(it)].permeability();
		}
		M_permeability /= M_aperture;
	}
}

#endif
