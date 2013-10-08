/*!
 *  @file interCellIntersections.hpp
 *  @brief Class to store intersections between a cell and a fault.
 *
 *  @author Luca Turconi <lturconi@gmail.com>
 *  @date 31-08-2012
 *
 */

#ifndef INTERCELLINTERSECTIONS_HPP_
#define INTERCELLINTERSECTIONS_HPP_

#include <utility>
#include <map>

#include "TypeDefinition.hpp"
#include "geomPoint3D.hpp"
#include "geomCPgridElements.hpp"
#include "geomCPgrid.hpp"

namespace Intersect
{
  /*!
    @typedef Intersection
      An intersection is defined as a pair.
      It stores:
       - the id of the edge on which intersection lies (UInt)
       - the intersection Point (Point3D)
    */
  typedef std::pair<UInt, Geometry::Point3D> Intersection;

  /*!
    @typedef CellIntersections_Iterator_Type
    Define iterator type on a CellIntersections object.
  */
  typedef std::map<UInt, Geometry::Point3D>::iterator      CellIntersections_Iterator_Type;

  /*!
    @typedef CellIntersections_Const_Iterator_Type
    Define const iterator type on a CellIntersections object.
  */
  typedef std::map<UInt, Geometry::Point3D>::const_iterator  CellIntersections_Const_Iterator_Type;

  /*!
    @class CellIntersections

    @author Luca Turconi <lturconi@gmail.com>

      This class defines an object to store cell intersections.

      The variables M_i, M_j, M_k identify the cell.
      The map M_intersections contains the intersections.

      Tools to compare two CellIntersections object and compute errors are available.

      Computed errors are stored in the M_intersectionsError map.
      It also provides some global statistics about errors.

      The method showMe allows to print the class informations.

    */
  class  CellIntersections
  {
  public:
    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    CellIntersections();

    //! Constructor, getting Corner Point cell
    /*!
     * @param c The Corner Point cell
     */
    CellIntersections (const Geometry::CPcell& c);

    //! Constructor, getting cell indices
    /*!
     * @param i The index i
     * @param j The index j
     * @param k The index k
     */
    CellIntersections (const UInt& i, const UInt& j, const UInt& k);

    //!Destructor
    ~CellIntersections();

    //@}

    //! @name Get Methods
    //@{

    //! Get index i (const)
    /*!
     * @return The index i
     */
    inline UInt i() const
    {
      return M_i;
    }

    //! Get index i
    /*!
     * @return A reference to the index i
     */
    inline UInt& i()
    {
      return M_i;
    }

    //! Get index j (const)
    /*!
     * @return The index j
     */
    inline UInt j() const
    {
      return M_j;
    }

    //! Get index j
    /*!
     * @return A reference to the index j
     */
    inline UInt& j()
    {
      return M_j;
    }

    //! Get index k (const)
    /*!
     * @return The index k
     */
    inline UInt k() const
    {
      return M_k;
    }

    //! Get index k
    /*!
     * @return A reference to the index k
     */
    inline UInt& k()
    {
      return M_k;
    }

    //! Get M_intersectionsError map
    /*!
     * @return A constant reference to the M_intersectionsError map
     */
    inline const std::map<UInt, Geometry::Point3D>& intersectionsErrorMap() const
    {
      return M_intersectionsError;
    }

    //! Get Mean Error
    /*!
     * @return The mean error
     */
    inline Real meanError() const
    {
      return M_meanError;
    }

    //! Get Max Error
    /*!
     * @return The max error
     */
    inline std::pair<UInt, Real> maxError() const
    {
      return M_maxError;
    }

    //! Get Missed Intersections counter
    /*!
     * Number of intersections stored only in exact solution.
     * (The current object doesn't store these intersections)
     * @return The number of missed intersections
     */
    inline UInt missedIntersections() const
    {
      return M_missedIntersections;
    }

    //! Get Additional Intersections counter
    /*!
     * Number of intersections stored only in the current object.
     * (The exact solution doesn't store these intersections)
     * @return The number of Additional intersections
     */
    inline UInt additionalIntersections() const
    {
      return M_additionalIntersections;
    }

    //@}

    //! @name Methods
    //@{

    //! Insert a new intersection
    /*!
     * @param inter The intersection to be inserted
     * @return TRUE -> operation ended correctly
          FALSE -> an error occurred
     */
    bool insert (const Intersection& inter, bool printError = 1);

    //! Clear the intersections map
    /*!
     * It removes all the elements from M_intersections map.
     */
    inline void clearIntersections()
    {
      M_intersections.clear();
    }

    //! Begin iterator
    /*!
     * @return The iterator to the first element of the M_intersections map.
     */
    inline CellIntersections_Iterator_Type begin()
    {
      return this->M_intersections.begin();
    }

    //! Begin iterator (const)
    /*!
     * @return The const iterator to the first element of the M_intersections map.
     */
    inline CellIntersections_Const_Iterator_Type begin() const
    {
      return this->M_intersections.begin();
    }

    //! End iterator
    /*!
     * @return The iterator to the last element of the M_intersections map.
     */
    inline CellIntersections_Iterator_Type end()
    {
      return this->M_intersections.end();
    }

    //! End iterator (const)
    /*!
     * @return The const iterator to the last element of the M_intersections map.
     */
    inline CellIntersections_Const_Iterator_Type end() const
    {
      return this->M_intersections.end();
    }

    //! M_intersections map size
    /*!
     * @return The size of M_intersections
     */
    inline UInt size() const
    {
      return M_intersections.size();
    }

    //! Find element
    /*!
     * Searches the container M_intersections for an element with id as key
     * returns an iterator to it if found, otherwise it returns an iterator to map::end
     * (the element past the end of the container).
     * @param id Key to be searched for
     * @return An iterator to the element, if the specified key value is found,
     *      or map::end if the specified key is not found in the container.
     */
    inline CellIntersections_Iterator_Type find ( const UInt& id )
    {
      return M_intersections.find (id);
    }

    //! Find element (const)
    /*!
     * Searches the container M_intersections for an element with id as key
     * returns a const iterator to it if found, otherwise it returns an iterator to map::end
     * (the element past the end of the container).
     * @param id Key to be searched for
     * @return A const iterator to the element, if the specified key value is found,
     *      or map::end if the specified key is not found in the container.
     */
    inline CellIntersections_Const_Iterator_Type find ( const UInt& id ) const
    {
      return M_intersections.find (id);
    }

    //! Import the intersections of cellInter in the current object
    /*!
     * @param cellInter The CellIntersections object to be merged
     * @return TRUE -> operation ended correctly
          FALSE -> an error occurred
     */
    bool importIntersections (const CellIntersections& cellInter);

    //! Compute errors
    /*!
     * It computes errors give an object of the same type storing the exact intersections.
     * Errors are saved in M_intersectionsError map.
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
    void computeErrors (const CellIntersections& exactInter);

    //! Display information in M_intersectionsError
    /*!
     * @param out Specify the output format (std::cout by default)
     */
    void showErrors (std::ostream&   out = std::cout) const;

    //! Display an error summary
    /*!
     * @param out Specify the output format (std::cout by default)
     */
    void errorSummary (std::ostream&   out = std::cout) const;

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
    UInt M_i, M_j, M_k;
    std::map<UInt, Geometry::Point3D> M_intersections;
    // Error tools
    std::map<UInt, Geometry::Point3D> M_intersectionsError;
    // (-1,NaN,NaN) -> intersection not found
    // (+1,NaN,NaN) -> found wrong intersection
    Real M_meanError;
    std::pair<UInt, Real> M_maxError;
    UInt M_missedIntersections;
    UInt M_additionalIntersections;
  };

} // namespace Intersect

#endif /* INTERCELLINTERSECTIONS_HPP_ */