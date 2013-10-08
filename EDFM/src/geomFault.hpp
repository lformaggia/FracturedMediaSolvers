/*!
*  @file geomFault.hpp
*  @brief Class for Fault.
*
*  @author Luca Turconi <lturconi@gmail.com>
*  @date 22-09-2012
*
*/

#ifndef GEOMFAULT_HPP_
#define GEOMFAULT_HPP_

#include "TypeDefinition.hpp"
#include "geomPoint3D.hpp"
#include "geomBilinearSurface.hpp"
#include "geomCPgridElements.hpp"
#include "geomCPgrid.hpp"
#include "interCellIntersections.hpp"
#include "interGridIntersections.hpp"
#include "interGridEdgeMap.hpp"
#include "interGridIntersectionMap.hpp"
#include "adtree.hpp"
#include "tolerances.hpp"
namespace Geometry
{
  /*!
    @class Fault

    @author Luca Turconi <lturconi@gmail.com>, Anna Scotti

      This class implements the concept of Fault.

      A fault is essentially a bilinear surface.
      It provides additional methods to compute intersections with Corner Point
      cells and Corner Point grids.

    */
  class Fault: public BilinearSurface
  {
  public:
    //! @name Constructor & Destructor
    //@{

    //! Empty constructor
    Fault();

    //! Constructor, getting the extremal points
    /*!
     * Attention: when creating the fault, a and c are considered as opposite points.
     * @param a The first point
     * @param b The second point
     * @param c The third point
     * @param d The fourth point
     */
    Fault (const Point3D& a, const Point3D& b, const Point3D& c, const Point3D& d);

    //! Copy constructor
    /*!
     * @param b The bilinear surface copied in the new object
     */
    Fault (const BilinearSurface& b);

    //! Destructor
    virtual ~Fault();



    //@}

    //! @name Methods
    //@{

    //! Approximated test for intersection with a cell
    /*!
     * The fault is approximated with two triangles.
     * The test for intersection is performed between the edges of the cell and this triangles.
     * @param c The cell to be tested
     * @param stdDivision A bool parameter to chose the way to split the fault.
     * @return TRUE if the cell has an intersection with a triangle
     *    FALSE if the intersection does not exist
     */
    bool approxIsIntersectedByCell (const CPcell& c, const bool& stdDivision = 1) const;

    //! Compute exact intersection with a cell
    /*!
     * It computes the intersection between the fault and the edges of a given cell.
     * The nonlinear problems are solved using the standard Newton's method.
     * If the intersection does not exist or in case of coplanar edges,
     * the object cellinter remains empty.
     * @param c The segment
     * @param cellinter The object to store intersections
     * @param toll The Newton's method tolerance
     * @param maxIter The Newton's method iteration limit
     */
    void newtonIntersectionWithCell
    (const CPcell& c, Intersect::CellIntersections& cellinter,
     const Real& toll = eps, const UInt& maxIter = 60) const;

    //! Compute approximated intersection with a cell
    /*!
     * The fault is approximated with two triangles.
     * It computes the intersection between the fault approximation and the edges of a given cell.
     * If the intersection does not exist or in case of coplanar edges,
     * the object cellinter remains empty.
     * @param c The segment
     * @param cellinter The object to store intersections
     * @param stdDivision A bool parameter to chose the way to split the fault.
     */
    void approxIntersectionWithCell
    (const CPcell& c, Intersect::CellIntersections& cellinter,
     const bool& stdDivision = 1) const;

    //! Approximated test with exact computation of intersection with a cell
    /*!
     * The fault is approximated with two triangles.
     * The method performs a test for the existence of intersections with the edges
     * of the cell using the surface approximation.
     * If the test return TRUE, it will compute the exact intersection between the fault
     * and the edges of a given cell.
     * The nonlinear problems are solved using the standard Newton's method.
     * If the intersection does not exist or in case of coplanar edges,
     * the object cellinter remains empty.
     * @param c The segment
     * @param cellinter The object to store intersections
     * @param toll The Newton's method tolerance
     * @param maxIter The Newton's method iteration limit
     */
    void approxNewtonIntersectionWithCell
    (const CPcell& c, Intersect::CellIntersections& cellinter,
     const Real& toll = eps, const UInt& maxIter = 60) const;

    //! Compute intersection between a given cell and the fault triangulation
    /*!
     * It computes the intersection between the fault triangulation
     * and the edges of a given cell.
     * If the intersection does not exist or in case of coplanar edges,
     * the object cellinter remains empty.
     * @param c The segment
     * @param cellinter The object to store intersections
     * @param toll The tolerance of the method
     */
    void approxRefinedIntersectionWithCell
    (const CPcell& c, Intersect::CellIntersections& cellinter,
     const Real& toll = eps);

    //! Compute exact intersection with a Corner Point grid
    /*!
     * It computes the intersection between the fault and a Corner Point grid.
     * The method uses 3 for loops to move through the cells of the grid.
     * The nonlinear problems are solved using the standard Newton's method.
     * If the intersection does not exist or in case of coplanar edges,
     * the object gridInter remains empty.
     * @param g The Corner Point grid
     * @param gridInter The object to store intersections
     * @param toll The Newton's method tolerance
     * @param maxIter The Newton's method iteration limit
     */
    void newtonIntersectionWithGrid_FOR3
    (const CPgrid& g, Intersect::GridIntersections& gridInter,
     const Real& toll = eps, const UInt& maxIter = 60) const;

    void newtonIntersectionWithGrid_FOR3OPT
    (const CPgrid& g, Intersect::GridIntersections& gridInter,
     const Real& toll = eps, const UInt& maxIter = 60) const;
    //! Computes approximated intersection with a Corner Point grid
    /*!
      It uses a bisection algorithm. It employs the new isIn2() method
      of CPCell, based on inverting the trilinear map.
     */

    void newtonIntersectionWithGrid_BISECTION
    (const CPgrid& g, Intersect::GridIntersections& gridInter,
     const Real& toll = eps, const UInt& maxIter = 60) const;

    //! Compute approximated intersection with a Corner Point grid
    /*!
     * The fault is approximated with two triangles.
     * It computes the intersection between the fault approximation and a Corner Point grid.
     * The method uses 3 for loops to move through the cells of the grid.
     * If the intersection does not exist or in case of coplanar edges,
     * the object cellinter remains empty.
     * @param g The Corner Point grid
     * @param gridInter The object to store intersections
     * @param stdDivision A bool parameter to chose the way to split the fault.
     */
    void approxIntersectionWithGrid_FOR3
    (const CPgrid& g, Intersect::GridIntersections& gridInter,
     const bool& stdDivision = 1) const;

    //! Approximated test with exact computation of intersection with a Corner Point grid
    /*!
     * The fault is approximated with two triangles.
     * The method performs a test for the existence of intersections with the
     * Corner Point grid using the surface approximation.
     * If the test return TRUE, it will compute the exact intersection between the fault
     * and the edges of a given cell.
     * The nonlinear problems are solved using the standard Newton's method.
     * If the intersection does not exist or in case of coplanar edges,
     * the object cellinter remains empty.
     * @param g The Corner Point grid
     * @param gridInter The object to store intersections
     * @param toll The Newton's method tolerance
     * @param maxIter The Newton's method iteration limit
     */
    void approxNewtonIntersectionWithGrid_FOR3
    (const CPgrid& g, Intersect::GridIntersections& gridInter,
     const Real& toll = eps, const UInt& maxIter = 60) const;

    //! Compute intersection between a Corner Point grid and the fault triangulation
    /*!
     * It computes the intersection between the fault triangulation
     * and a Corner Point grid.
     * If the intersection does not exist or in case of coplanar edges,
     * the object cellinter remains empty.
     * @param g The Corner Point grid
     * @param gridInter The object to store intersections
     * @param toll The tolerance of the method
     */
    void approxRefinedIntersectionWithGrid_FOR3
    (const CPgrid& g, Intersect::GridIntersections& gridInter,
     const Real& toll = eps);

    //! Compute exact intersection with a GridEdgeMap object
    /*!
     * It computes the intersection between the fault and the edges of a Corner Point grid.
     * The method uses the GridEdgeMap class to store grid edges without repetitions.
     * The nonlinear problems are solved using the standard Newton's method.
     * If the intersection does not exist or in case of coplanar edges,
     * the object gridInter remains empty.
     * @param gridEdge The Corner Point grid edges (without repetitions)
     * @param gridInter The object to store intersections (without repetitions)
     * @param toll The Newton's method tolerance
     * @param maxIter The Newton's method iteration limit
     */
    void newtonIntersectionWithGridEdge
    (const Intersect::GridEdgeMap& gridEdge, Intersect::GridIntersectionMap& gridInter,
     const Real& toll = eps, const UInt& maxIter = 60) const;

    //! Compute approximated intersection with a GridEdgeMap object
    /*!
     * The fault is approximated with two triangles.
     * It computes the intersection between the fault approximation and the edges of a Corner Point grid.
     * The method uses the GridEdgeMap class to store grid edges without repetitions.
     * If the intersection does not exist or in case of coplanar edges,
     * the object cellinter remains empty.
     * @param gridEdge The Corner Point grid edges (without repetitions)
     * @param gridInter The object to store intersections (without repetitions)
     * @param stdDivision A bool parameter to chose the way to split the fault.
     */
    void approxIntersectionWithGridEdge
    (const Intersect::GridEdgeMap& gridEdge, Intersect::GridIntersectionMap& gridInter,
     const bool& stdDivision = 1) const;

    //! Approximated test with exact computation of intersection with a GridEdgeMap object
    /*!
     * The fault is approximated with two triangles.
     * The method performs a test for the existence of intersections with the edges of a
     * Corner Point grid using the surface approximation.
     * If the test return TRUE, it will compute the exact intersection.
     * The nonlinear problems are solved using the standard Newton's method.
     * The method uses the GridEdgeMap class to store grid edges without repetitions.
     * If the intersection does not exist or in case of coplanar edges,
     * the object cellinter remains empty.
     * @param gridEdge The Corner Point grid edges (without repetitions)
     * @param gridInter The object to store intersections (without repetitions)
     * @param toll The Newton's method tolerance
     * @param maxIter The Newton's method iteration limit
     */
    void approxNewtonIntersectionWithGridEdge
    (const Intersect::GridEdgeMap& gridEdge, Intersect::GridIntersectionMap& gridInter,
     const Real& toll = eps, const UInt& maxIter = 60) const;

    //! Compute intersection between a GridEdgeMap object and the fault triangulation
    /*!
     * It computes the intersection between the fault triangulation
     * and the edges of a Corner Point grid.
     * The method uses the GridEdgeMap class to store grid edges without repetitions.
     * If the intersection does not exist or in case of coplanar edges,
     * the object cellinter remains empty.
     * @param gridEdge The Corner Point grid edges (without repetitions)
     * @param gridInter The object to store intersections (without repetitions)
     * @param toll The tolerance of the method
     */
    void approxRefinedIntersectionWithGridEdge
    (const Intersect::GridEdgeMap& gridEdge, Intersect::GridIntersectionMap& gridInter,
     const Real& toll = eps);

    //! Compute exact intersection with a Corner Point grid
    /*!
     * It computes the intersection between the fault and a Corner Point grid.
     * The method uses the GridEdgeMap class to store grid edges without repetitions.
     * It is essentially based on newtonIntersectionWithGridEdge() method.
     * The nonlinear problems are solved using the standard Newton's method.
     * If the intersection does not exist or in case of coplanar edges,
     * the object gridInter remains empty.
     * @param g The Corner Point grid
     * @param gridInter The object to store intersections
     * @param toll The Newton's method tolerance
     * @param maxIter The Newton's method iteration limit
     */
    void newtonIntersectionWithGrid_EDGEMAP
    (const CPgrid& g, Intersect::GridIntersections& gridInter,
     const Real& toll = eps, const UInt& maxIter = 60) const;

    //! Compute approximated intersection with a Corner Point grid
    /*!
     * The fault is approximated with two triangles.
     * It computes the intersection between the fault approximation and a Corner Point grid.
     * The method uses the GridEdgeMap class to store grid edges without repetitions.
     * It is essentially based on approxIntersectionWithGridEdge() method.
     * If the intersection does not exist or in case of coplanar edges,
     * the object cellinter remains empty.
     * @param g The Corner Point grid
     * @param gridInter The object to store intersections
     * @param stdDivision A bool parameter to chose the way to split the fault.
     */
    void approxIntersectionWithGrid_EDGEMAP
    (const CPgrid& g, Intersect::GridIntersections& gridInter,
     const bool& stdDivision = 1) const;

    //! Approximated test with exact computation of intersection with a Corner Point grid
    /*!
     * The fault is approximated with two triangles.
     * The method performs a test for the existence of intersections with the
     * Corner Point grid using the surface approximation.
     * If the test return TRUE, it will compute the exact intersection between the fault
     * and the edges of a given cell.
     * The nonlinear problems are solved using the standard Newton's method.
     * It is essentially based on approxNewtonIntersectionWithGridEdge() method.
     * If the intersection does not exist or in case of coplanar edges,
     * the object cellinter remains empty.
     * @param g The Corner Point grid
     * @param gridInter The object to store intersections
     * @param toll The Newton's method tolerance
     * @param maxIter The Newton's method iteration limit
     */
    void approxNewtonIntersectionWithGrid_EDGEMAP
    (const CPgrid& g, Intersect::GridIntersections& gridInter,
     const Real& toll = eps, const UInt& maxIter = 60) const;

    //! Compute intersection between a Corner Point grid and the fault triangulation
    /*!
     * It computes the intersection between the fault triangulation
     * and a Corner Point grid.
     * The method uses the GridEdgeMap class to store grid edges without repetitions.
     * It is essentially based on approxRefinedIntersectionWithGridEdge() method.
     * If the intersection does not exist or in case of coplanar edges,
     * the object cellinter remains empty.
     * @param g The Corner Point grid
     * @param gridInter The object to store intersections
     * @param toll The tolerance of the method
     */
    void approxRefinedIntersectionWithGrid_EDGEMAP
    (const CPgrid& g, Intersect::GridIntersections& gridInter,
     const Real& toll = eps);

    //@}
  };

} // namespace Geometry

#endif /* GEOMFAULT_HPP_ */
