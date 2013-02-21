 /*!
 *	@file interGridIntersections.hpp
 *	@brief Class to store intersections between a grid and faults.
 *
 *	@author Luca Turconi <lturconi@gmail.com>
 *  @date 31-08-2012
 *
 */

 
#ifndef INTERGRIDINTERSECTIONS_HPP_
#define INTERGRIDINTERSECTIONS_HPP_

#include <utility>
#include <map>

#include "TypeDefinition.hpp"
#include "geomPoint3D.hpp"
#include "geomCPgridElements.hpp"
#include "geomCPgrid.hpp"
#include "interCellIntersections.hpp"

namespace Intersect
{
	/*!
		@typedef GridIntersections_Iterator_Type
		Define iterator type on a GridIntersections object.
	*/
typedef std::map<UInt,CellIntersections>::iterator			GridIntersections_Iterator_Type;

	/*!
		@typedef GridIntersections_Const_Iterator_Type
		Define const iterator type on a GridIntersections object.
	*/
typedef std::map<UInt,CellIntersections>::const_iterator	GridIntersections_Const_Iterator_Type;

	/*!
		@class GridIntersections
		
		@author Luca Turconi <lturconi@gmail.com>
		
    	This class defines an object to store grid intersections.
    	
    	The variables M_Nx, M_Ny, M_Nz set the grid dimensions.
    	Intersections data are stored in M_CellIntersectionsMap.
    	The map key follows the rule:
			ID = i + (j-1)*M_Nx + (k-1)*M_Nx*M_Ny
    	
    	Tools to compare two GridIntersections object and compute errors are available.
    	It also provides some global statistics about errors.
    	
    	The method showMe allows to print the class informations.

    */
class GridIntersections
{
public:
	//! @name Constructors & Destructor
	//@{
	
	//! Constructor, getting a Corner Point grid
	/*!
	 * @param g The Corner Point grid
	 */
	GridIntersections(const Geometry::CPgrid & g);
	
	//! Constructor, getting the grid dimensions
	/*!
	 * @param Nx The number of cells in x direction
	 * @param Ny The number of cells in y direction
	 * @param Nz The number of cells in z direction
	 */
	GridIntersections(const UInt & Nx, const UInt & Ny, const UInt & Nz);
	
	//! Destructor
	~GridIntersections();
	
	//@}
	
	//! @name Get Methods
	//@{
		
	//! Get x dimension
	/*!
	 * @return The number of cells in x direction
	 */
	inline UInt Nx() const { return M_Nx; }
	
	//! Get y dimension
	/*!
	 * @return The number of cells in y direction
	 */
	inline UInt Ny() const { return M_Ny; }
	
	//! Get z dimension
	/*!
	 * @return The number of cells in z direction
	 */
	inline UInt Nz() const { return M_Nz; }
	
	//! Get the intersections counter
	/*!
	 * @return The number of intersections
	 */
	inline UInt Nintersections() const { return M_Nintersections; }
	
	//! Get M_CellIntersectionsMap
	/*!
	 * @return A constant reference to the cell intersections map
	 */
	inline const std::map<UInt, CellIntersections> & cellIntersectionsMap() const
		{ return M_CellIntersectionsMap; }
	
	//! Get cell index i from map key
	/*!
	 * @param id The key value of the cell in the map
	 * @return The cell index i
	 */
	inline UInt i(const UInt & id) const
		{ return id!=0 ? (id-1)%M_Nx + 1 : 0; }
	
	//! Get cell index j from map key
	/*!
	 * @param id The key value of the cell in the map
	 * @return The cell index j
	 */
	inline UInt j(const UInt & id) const
		{ return id!=0 ? ((id-1)/M_Nx)%M_Ny +1 : 0; }
	
	//! Get cell index k from map key
	/*!
	 * @param id The key value of the cell in the map
	 * @return The cell index k
	 */
	inline UInt k(const UInt & id) const
		{ return id!=0 ? ((id-1)/M_Nx)/M_Ny +1 : 0; }
	
	//! Get CellIntersections object from map key
	/*!
	 * @param id The key value of the cell in the map
	 * @return The CellIntersections object
	 */
	CellIntersections getCellIntersections(const UInt & id) const;
	
	//! Get CellIntersections object from (i,j,k) indexes
	/*!
	 * @param i The index i
	 * @param j The index j
	 * @param k The index k
	 * @return The CellIntersections object
	 */
	CellIntersections getCellIntersections(const UInt & i, const UInt & j, const UInt & k) const;
	
	//! Get Mean Error
	/*!
	 * @return The mean error
	 */
	inline Real meanError() const
		{ return M_meanError; }
	
	//! Get Max Error
	/*!
	 * @return The max error
	 */
	inline std::pair<UInt,Real> maxError() const
		{ return M_maxError; }

	//! Get Missed Cells counter
	/*!
	 * Number of cells stored only in exact solution.
	 * (The current object doesn't store these cells)
	 * @return The number of missed cells
	 */
	inline UInt missedCells() const
		{ return M_missedCells; }
	
	//! Get Additional Cells counter
	/*!
	 * Number of cells stored only in the current object.
	 * (The exact solution doesn't store these cells)
	 * @return The number of additional cells
	 */
	inline UInt additionalCells() const
		{ return M_additionalCells; }
	
	//@}
	
	//! @name Set Methods
	//@{
	
	//! Set the intersections counter
	/*!
	 * @param 
	 */
	inline void setNintersections(const UInt & n)
		{ M_Nintersections = n; }
	
	//! Increment the intersections counter
	/*!
	 * @param 
	 */
	inline void incrementNintersections()
		{ M_Nintersections +=1; }
	
	//@}
	
	//! @name Methods
	//@{
	
	//! Insert a new CellIntersections object
	/*!
	 * @param inter The CellIntersections object to be insert
	 * @return TRUE -> operation ended correctly
				FALSE -> an error occurred
	 */
	bool insert(const CellIntersections & inter, bool printError=1);
	
	//! M_CellIntersectionsMap map size
	/*!
	 * @return The size of M_CellIntersectionsMap
	 */
	inline std::map<UInt,CellIntersections>::size_type size() const
		{ return M_CellIntersectionsMap.size(); }
	
	//! Begin iterator
	/*!
	 * @return The iterator to the first element of the M_CellIntersectionsMap.
	 */
	inline GridIntersections_Iterator_Type begin()
		{ return M_CellIntersectionsMap.begin(); }
		
	//! Begin iterator (const)
	/*!
	 * @return The const iterator to the first element of the M_CellIntersectionsMap.
	 */
	inline GridIntersections_Const_Iterator_Type begin() const
		{ return M_CellIntersectionsMap.begin(); }
	
	//! End iterator
	/*!
	 * @return The iterator to the last element of the M_CellIntersectionsMap.
	 */
	inline GridIntersections_Iterator_Type end()
		{ return M_CellIntersectionsMap.end(); }
		
	//! End iterator (const)
	/*!
	 * @return The const iterator to the last element of the M_CellIntersectionsMap.
	 */
	inline GridIntersections_Const_Iterator_Type end() const
		{ return M_CellIntersectionsMap.end(); }
	
	//! Find element
	/*!
	 * Searches the container M_CellIntersectionsMap for an element with id as key
	 * returns an iterator to it if found, otherwise it returns an iterator to map::end
	 * (the element past the end of the container).
	 * @param id Key to be searched for
	 * @return An iterator to the element, if the specified key value is found,
	 * 			or map::end if the specified key is not found in the container.
	 */
	inline GridIntersections_Iterator_Type find ( const UInt & id )
		{ return M_CellIntersectionsMap.find(id); }
	
	//! Find element (const)
	/*!
	 * Searches the container M_CellIntersectionsMap for an element with id as key
	 * returns a const iterator to it if found, otherwise it returns an iterator to map::end
	 * (the element past the end of the container).
	 * @param id Key to be searched for
	 * @return A const iterator to the element, if the specified key value is found,
	 * 			or map::end if the specified key is not found in the container.
	 */
	inline GridIntersections_Const_Iterator_Type find ( const UInt & id ) const
		{ return M_CellIntersectionsMap.find(id); }
	
	//! Erase element (iterator)
	/*!
	 * Removes from the container M_CellIntersectionsMap the element pointed by position.
	 * @param position The iterator pointing the element to be erased.
	 */
	inline void erase ( GridIntersections_Iterator_Type position )
		{ 
			M_Nintersections -= (*position).second.size();
			M_CellIntersectionsMap.erase(position);
		}
	
	//! Erase element (key)
	/*!
	 * Removes from the container M_CellIntersectionsMap the element with the key value of x.
	 * @param x The key of the element to be erased.
	 * @return The number of elements erased, which in map containers is 1
	 * 			if an element with a key value of x existed
	 * 			(and thus was subsequently erased), and zero otherwise.
	 */
	inline std::map<UInt,CellIntersections>::size_type erase ( const UInt& x )
		{ 
			M_Nintersections -= M_CellIntersectionsMap.find(x)->second.size();
			return M_CellIntersectionsMap.erase(x);
		}
		
	//! Erase elements (range)
	/*!
	 * Removes from the container M_CellIntersectionsMap a range of elements ([first,last)).
	 * @param first The iterator pointing the first element to be erased.
	 * @param last The iterator pointing the last element to be erased.
	 */
    inline void erase ( GridIntersections_Iterator_Type first, GridIntersections_Iterator_Type last )
		{ 
			for(GridIntersections_Iterator_Type it=first; it!=last; ++it)
				M_Nintersections -= (*it).second.size();
			M_CellIntersectionsMap.erase(first,last);
		}
	
	//! Clear the intersections map
	/*!
	 * It removes all the elements from M_CellIntersectionsMap.
	 */
	inline void clearMap()
		{ 
			M_Nintersections=0;
			M_CellIntersectionsMap.clear();
		}
	
	//! Clear GridIntersections
	/*!
	 * It clears all data in GridIntersections.
	 * Only the variables M_Nx, M_Ny, M_Nz are preserved.
	 */
	void clearAll();
	
	//! Compute errors
	/*!
	 * It computes errors give an object of the same type storing the exact intersections.
	 * 
	 * For points stored in the current object but not in the exact solution
	 * (found wrong intersection)
	 * the error map contains the element (+1,NaN,NaN).
	 * 
	 * For points stored in the exact solution but not in the current object
	 * (intersection not found):
	 * the error map contains the element (-1,NaN,NaN).
	 * 
	 * The method generates the global statistics about errors.
	 * @param exactInter The object storing the exact intersections.
	 */
	void computeErrors(const GridIntersections & exactInter);
	
	//! Display an error summary
	/*!
	 * @param out Specify the output format (std::cout by default)
	 */
	void errorSummary(std::ostream  & out=std::cout) const;
	
	//! Export in vtk format
	/*!
	 * It generates the vtk file, for the 3D visualization of the intersections.
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
	
	//@}
	
private:
	std::map<UInt, CellIntersections> M_CellIntersectionsMap;
	// ID = i + (j-1)*M_Nx + (k-1)*M_Nx*M_Ny
	UInt M_Nx, M_Ny, M_Nz;
	UInt M_Nintersections;
	// Error tools
	Real M_meanError;
	std::pair<UInt,Real> M_maxError;
	UInt M_missedCells;
	UInt M_additionalCells;
};

} // namespace Intersect

#endif /* INTERGRIDINTERSECTIONS_HPP_ */
