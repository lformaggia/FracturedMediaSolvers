/*!
 *	@file interGridEdgeMap.hpp
 *	@brief Class for storage of grid edges without repetitions.
 *
 *	@author Luca Turconi <lturconi@gmail.com>
 *  @date 27-09-2012
 *
 */


#ifndef INTERGRIDEDGEMAP_HPP_
#define INTERGRIDEDGEMAP_HPP_

#include <iostream>
#include <utility>
#include <map>
#include <list>
#include <cmath>
#include <functional>

#include "TypeDefinition.hpp"
#include "geomPoint3D.hpp"
#include "geomSegment.hpp"
#include "geomCPgrid.hpp"


namespace Intersect
{
	/*!
		@typedef LinkerElement
    	Define the type for the Linker elements.
    */
typedef std::pair<UInt,UInt> LinkerElement;

	/*!
		@typedef Linker_Iterator_Type
    	Define iterator type on a Linker object.
    */
typedef std::list<LinkerElement>::iterator Linker_Iterator_Type;

	/*!
		@typedef Linker_Const_Iterator_Type
    	Define const iterator type on a Linker object.
    */
typedef std::list<LinkerElement>::const_iterator Linker_Const_Iterator_Type;

	/*!
		@class Linker
		
		@author Luca Turconi <lturconi@gmail.com>
		
    	A Linker is a map storing unsigned int pairs. \n
    	It stores the connections between an edge
    	and the positions it occupies in a Corner Point grid. \n
		The map key rapresents the ID of the cells from which the edge comes
		(ID = i + (j-1)*M_Nx + (k-1)*M_Nx*M_Ny). \n
    	The associated UInt is the edge ID (from 1 to 12)
    	which identifies its internal position in the cell.

    */
class Linker
{
public:
	//! @name Constructors & Destructor
	//@{
		
	//! Empty Constructor
	Linker();
	
	//! Copy Constructor
	/*!
	 * @param l The Linker copied in the new object
	 */
	Linker(const Linker & l);
	
	//!Destructor
	~Linker();
	
	//@}
	
	//! @name Get Methods
	//@{
	
	//! Get M_LinksList reference
	/*!
	 * @return A constant reference to M_LinksList
	 */
	inline const std::list<LinkerElement> & getLinksList() const
		{ return M_LinksList; }
	
	//@}
	
	//! @name Methods
	//@{
	
	//! Insert a new link at the beginning of the list
	/*!
	 * @param link The link to be inserted
	 * @return TRUE -> operation ended correctly
				FALSE -> an error occurred
	 */
	inline void insert(const LinkerElement & link)
		{ M_LinksList.push_front(link); }
	
	//! Merge the element of linker l in the current linker
	/*!
	 * @param l The linker to be merged
	 * @return TRUE -> operation ended correctly
				FALSE -> an error occurred
	 */
	void merge(const Linker & l);
	
	//! Clear links list
	/*!
	 * It removes all the elements from M_LinksList.
	 */
	inline void clearLinksList()
		{ M_LinksList.clear(); }
	
	//! Begin iterator
	/*!
	 * @return The iterator to the first element of the M_LinksList.
	 */
	inline Linker_Iterator_Type begin()
		{ return this->M_LinksList.begin(); }
		
	//! Begin iterator (const)
	/*!
	 * @return The const iterator to the first element of the M_LinksList.
	 */
	inline Linker_Const_Iterator_Type begin() const
		{ return this->M_LinksList.begin(); }
	
	//! End iterator
	/*!
	 * @return The iterator to the last element of the M_LinksList.
	 */
	inline Linker_Iterator_Type end()
		{ return this->M_LinksList.end(); }
		
	//! End iterator (const)
	/*!
	 * @return The const iterator to the last element of the M_LinksList.
	 */
	inline Linker_Const_Iterator_Type end() const
		{ return this->M_LinksList.end(); }
	
	//! Display general information about the content of the class
	/*!
	 * List of things displayed in the class
	 * @param out Specify the output format (std::cout by default)
	 */
	void showMe(std::ostream  & out=std::cout) const;
	
	//@}
private:
	std::list<LinkerElement> M_LinksList;
};

	/*!
		@typedef GridEdgeMapElement
    	Define the value_type for the container GridEdgeMap.
    */
typedef std::pair<Geometry::Segment,Linker> GridEdgeMapElement;

	/*!
		@typedef GridEdgeMap_Iterator_Type
		Define iterator type on a GridEdgeMap object.
	*/
typedef std::map<Geometry::Segment,Linker>::iterator GridEdgeMap_Iterator_Type;

	/*!
		@typedef GridEdgeMap_Const_Iterator_Type
		Define const iterator type on a GridEdgeMap object.
	*/
typedef std::map<Geometry::Segment,Linker>::const_iterator GridEdgeMap_Const_Iterator_Type;

	/*!
		@class GridEdgeMap
		
		@author Luca Turconi <lturconi@gmail.com>
		
    	This class defines an object to store grid edges without repetitions.\n
    	\n
    	The variables M_Nx, M_Ny, M_Nz define the grid dimensions.\n
    	Edges are stored in M_EdgeMap.
    	Their positions in the Corner Point grid are described by the Linker.\n
    	\n
    	The method showMe allows to print the class informations.

    */
class GridEdgeMap
{
public:
	//! @name Constructors & Destructor
	//@{
		
	//! Empty Constructor
	GridEdgeMap();
	
	//! Constructor, getting the grid dimensions
	/*!
	 * @param Nx The number of cells in x direction
	 * @param Ny The number of cells in y direction
	 * @param Nz The number of cells in z direction
	 */
	GridEdgeMap(const UInt & Nx, const UInt & Ny, const UInt & Nz);
	
	//! Constructor, getting Corner Point grid
	/*!
	 * @param grid The Corner Point grid
	 */
	GridEdgeMap(const Geometry::CPgrid & grid);
	
	//! Copy Constructor
	/*!
	 * @param gem The GridEdgeMap copied in the new object
	 */
	GridEdgeMap(const GridEdgeMap & gem);
		
	//!Destructor
	~GridEdgeMap();
	
	//@}
	
	//! @name Get Methods
	//@{
		
	//! Get grid x dimension
	/*!
	 * @return The number of cells in x direction
	 */
	inline UInt Nx() const { return M_Nx; }
	
	//! Get grid y dimension
	/*!
	 * @return The number of cells in y direction
	 */
	inline UInt Ny() const { return M_Ny; }
	
	//! Get grid z dimension
	/*!
	 * @return The number of cells in z direction
	 */
	inline UInt Nz() const { return M_Nz; }
	
	//! Get M_EdgeMap reference
	/*!
	 * @return A constant reference to M_EdgeMap
	 */
	inline const std::map<Geometry::Segment,Linker> & getEdgeMap() const
		{ return M_EdgeMap; }
	
	//@}
	
	//! @name Methods
	//@{
		
	//! Extract edges from grid
	/*!
	 * Fill the current GridEdgeMap object with edges of the Corner Point grid g. \n
	 * Warning: dimensions of grid and GridEdgeMap must agree!
	 * @param g The Corner Point grid
	 * @return TRUE -> operation ended correctly
				FALSE -> an error occurred
	 */
	bool extractEdges(const Geometry::CPgrid & g);
		
	//! Insert a new edge
	/*!
	 * @param mapElem The Edge to be inserted
	 * @return TRUE -> operation ended correctly
				FALSE -> an error occurred
	 */
	void insert(const GridEdgeMapElement & mapElem);
	
	//! Clear the edge map
	/*!
	 * It removes all the elements from M_EdgeMap.
	 */
	inline void clearEdgeMap()
		{ M_EdgeMap.clear(); }
	
	//! Begin iterator
	/*!
	 * @return The iterator to the first element of the M_EdgeMap.
	 */
	inline GridEdgeMap_Iterator_Type begin()
		{ return this->M_EdgeMap.begin(); }
	
	//! Begin iterator (const)
	/*!
	 * @return The const iterator to the first element of the M_EdgeMap.
	 */
	inline GridEdgeMap_Const_Iterator_Type begin() const
		{ return this->M_EdgeMap.begin(); }
	
	//! End iterator
	/*!
	 * @return The iterator to the last element of the M_EdgeMap.
	 */
	inline GridEdgeMap_Iterator_Type end()
		{ return this->M_EdgeMap.end(); }
	
	//! End iterator (const)
	/*!
	 * @return The const iterator to the last element of the M_EdgeMap.
	 */
	inline GridEdgeMap_Const_Iterator_Type end() const
		{ return this->M_EdgeMap.end(); }
	
	//! M_EdgeMap map size
	/*!
	 * @return The size of M_EdgeMap
	 */
	inline UInt size() const
		{ return M_EdgeMap.size(); }
	
	//! Find element
	/*!
	 * Searches the container M_EdgeMap for a given Segment and
	 * returns an iterator to it if found, otherwise it returns an iterator to map::end
	 * (the element past the end of the container).
	 * @param s Segment to be searched for
	 * @return An iterator to the element, if the specified key value is found,
	 * 			or map::end if the specified key is not found in the container.
	 */
	inline GridEdgeMap_Iterator_Type find ( const Geometry::Segment & s )
		{ return M_EdgeMap.find(s); }
	
	//! Find element (const)
	/*!
	 * Searches the container M_EdgeMap for a given Segment and
	 * returns an iterator to it if found, otherwise it returns an iterator to map::end
	 * (the element past the end of the container).
	 * @param s Segment to be searched for
	 * @return A const iterator to the element, if the specified key value is found,
	 * 			or map::end if the specified key is not found in the container.
	 */
	inline GridEdgeMap_Const_Iterator_Type find ( const Geometry::Segment & s ) const
		{ return M_EdgeMap.find(s); }
	
	//! Export in vtk format
	/*!
	 * It generates the vtk file, for the 3D visualization of the edges.
	 * (Use paraview to open the vtk file)
	 * @param filename The name of the vtk file created by this method
	 * @return TRUE -> operation ended correctly
				FALSE -> an error occurred
	 */	
	bool exportVtk(const std::string & filename) const;
	
	//! Display general information about the content of the class
	/*!
	 * List of things displayed in the class
	 * @param out Specify the output format (std::cout by default)
	 */
	void showMe(std::ostream  & out=std::cout) const;
	
private:
	UInt M_Nx, M_Ny, M_Nz;
	std::map<Geometry::Segment,Linker> M_EdgeMap;
};

	/*!
		@fn Iterator safe_advancer
    	This function allows to advance a generic iterator of delta steps.
    	The safeness means that it always stops the increment when the end is reached.
    	@param it The iterator to be advanced
    	@param end The end limit for advancing operation
    	@param delta The steps of single advance operation
    	@return The advanced iterator
    */
template <typename Itr>
Itr safe_advancer(Itr it, Itr end, UInt delta)
{
    while(it != end && delta--)
        it++;
    return it;
}

} // namespace Intersect

#endif /* INTERGRIDEDGEMAP_HPP_ */