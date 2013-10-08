/*!
*  @file interGridIntersectionMap.hpp
*  @brief Class for storage of grid intersections without repetitions.
*
*  @author Luca Turconi <lturconi@gmail.com>
*  @date 27-09-2012
*
*/


#ifndef INTERGRIDINTERSECTIONMAP_HPP
#define INTERGRIDINTERSECTIONMAP_HPP

#include <iostream>
#include <utility>
#include <map>

#include "TypeDefinition.hpp"
#include "geomPoint3D.hpp"
#include "geomCPgrid.hpp"
#include "interGridEdgeMap.hpp"
#include "interGridIntersections.hpp"

namespace Intersect
{
  /*!
    @typedef GridIntersectionMapElement
      Define the value_type for the container GridIntersectionMapElement.
    */
  typedef std::pair<Geometry::Point3D, Linker> GridIntersectionMapElement;

  /*!
    @typedef GridIntersectionMap_Iterator_Type
    Define iterator type on a GridIntersectionMap object.
  */
  typedef std::map<Geometry::Point3D, Linker>::iterator GridIntersectionMap_Iterator_Type;

  /*!
    @typedef GridIntersectionMap_Const_Iterator_Type
    Define const iterator type on a GridIntersectionMap object.
  */
  typedef std::map<Geometry::Point3D, Linker>::const_iterator GridIntersectionMap_Const_Iterator_Type;

  /*!
    @class GridIntersectionMap

    @author Luca Turconi <lturconi@gmail.com>

      This class defines an object to store grid intersections without repetitions.\n
      \n
      The variables M_Nx, M_Ny, M_Nz define the grid dimensions.\n
      Intersections are stored in M_IntersectionMap.
      Their positions in the Corner Point grid are described by the Linker.\n
      \n
      The method showMe allows to print the class informations.

    */
  class GridIntersectionMap
  {
  public:
    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    GridIntersectionMap();

    //! Constructor, getting the grid dimensions
    /*!
     * @param Nx The number of cells in x direction
     * @param Ny The number of cells in y direction
     * @param Nz The number of cells in z direction
     */
    GridIntersectionMap (const UInt& Nx, const UInt& Ny, const UInt& Nz);

    //! Constructor, getting Corner Point grid
    /*!
     * @param grid The Corner Point grid
     */
    GridIntersectionMap (const Geometry::CPgrid& grid);

    //! Copy Constructor
    /*!
     * @param gem The GridIntersectionMap copied in the new object
     */
    GridIntersectionMap (const GridIntersectionMap& gim);

    //!Destructor
    ~GridIntersectionMap();

    //@}

    //! @name Get Methods
    //@{

    //! Get grid x dimension
    /*!
     * @return The number of cells in x direction
     */
    inline UInt Nx() const
    {
      return M_Nx;
    }

    //! Get grid y dimension
    /*!
     * @return The number of cells in y direction
     */
    inline UInt Ny() const
    {
      return M_Ny;
    }

    //! Get grid z dimension
    /*!
     * @return The number of cells in z direction
     */
    inline UInt Nz() const
    {
      return M_Nz;
    }

    //! Get M_IntersectionMap reference
    /*!
     * @return A constant reference to M_IntersectionMap
     */
    inline const std::map<Geometry::Point3D, Linker>& getIntersectionMap() const
    {
      return M_IntersectionMap;
    }

    //@}

    //! @name Methods
    //@{

    //! Insert a new intersection
    /*!
     * @param mapElem The intersection to be inserted
     * @return TRUE -> operation ended correctly
          FALSE -> an error occurred
     */
    void insert (const GridIntersectionMapElement& mapElem);

    //! Clear the intersections map
    /*!
     * It removes all the elements from M_IntersectionMap.
     */
    inline void clearIntersectionMap()
    {
      M_IntersectionMap.clear();
    }

    //! Begin iterator
    /*!
     * @return The iterator to the first element of the M_IntersectionMap.
     */
    inline GridIntersectionMap_Iterator_Type begin()
    {
      return this->M_IntersectionMap.begin();
    }

    //! Begin iterator (const)
    /*!
     * @return The const iterator to the first element of the M_IntersectionMap.
     */
    inline GridIntersectionMap_Const_Iterator_Type begin() const
    {
      return this->M_IntersectionMap.begin();
    }

    //! End iterator
    /*!
     * @return The iterator to the last element of the M_IntersectionMap.
     */
    inline GridIntersectionMap_Iterator_Type end()
    {
      return this->M_IntersectionMap.end();
    }

    //! End iterator (const)
    /*!
     * @return The const iterator to the last element of the M_IntersectionMap.
     */
    inline GridIntersectionMap_Const_Iterator_Type end() const
    {
      return this->M_IntersectionMap.end();
    }

    //! M_IntersectionMap map size
    /*!
     * @return The size of M_IntersectionMap
     */
    inline UInt size() const
    {
      return M_IntersectionMap.size();
    }

    //! Find element
    /*!
     * Searches the container M_IntersectionMap for a given Point and
     * returns an iterator to it if found, otherwise it returns an iterator to map::end
     * (the element past the end of the container).
     * @param p Point to be searched for
     * @return An iterator to the element, if the specified key value is found,
     *      or map::end if the specified key is not found in the container.
     */
    inline GridIntersectionMap_Iterator_Type find ( const Geometry::Point3D& p )
    {
      return M_IntersectionMap.find (p);
    }

    //! Find element (const)
    /*!
     * Searches the container M_IntersectionMap for a given Point and
     * returns an iterator to it if found, otherwise it returns an iterator to map::end
     * (the element past the end of the container).
     * @param p Point to be searched for
     * @return A const iterator to the element, if the specified key value is found,
     *      or map::end if the specified key is not found in the container.
     */
    inline GridIntersectionMap_Const_Iterator_Type find ( const Geometry::Point3D& p ) const
    {
      return M_IntersectionMap.find (p);
    }

    //! Export data from GridIntersectionMap to GridIntersections
    /*!
     * @param gridInter The target GridIntersections object
     * @return TRUE -> operation ended correctly
          FALSE -> an error occurred
     */
    bool exportToGridIntersections (GridIntersections& gridInter) const;

    //! Import data from GridIntersections to GridIntersectionMap
    /*!
     * @param gridInter The target GridIntersections object
     * @return TRUE -> operation ended correctly
          FALSE -> an error occurred
     */

    bool importFromGridIntersections (GridIntersections& gridInter);

    //! Export in vtk format
    /*!
     * It generates the vtk file, for the 3D visualization of the intersections.
     * (Use paraview to open the vtk file)
     * @param filename The name of the vtk file created by this method
     * @return TRUE -> operation ended correctly
          FALSE -> an error occurred
     */
    bool exportVtk (const std::string& filename) const;

    //! Display general information about the content of the class
    /*!
     * List of things displayed in the class
     * @param out Specify the output format (std::cout by default)
     */
    void showMe (std::ostream&   out = std::cout) const;

  private:
    UInt M_Nx, M_Ny, M_Nz;
    std::map<Geometry::Point3D, Linker> M_IntersectionMap;
  };

} // namespace Intersect

#endif /* INTERGRIDINTERSECTIONMAP_HPP */