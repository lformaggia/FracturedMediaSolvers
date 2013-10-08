/*!
 *  @file cutCellProperies.hpp
 *  @author Luca Turconi <lturconi@gmail.com>
 *  @author Anna Scotti <annascotti@hotmail.com>
 *  @author Luca Formaggia <luca.formaggia@gmail.com>
 *  @date 22-08-2013
 *
 */

#ifndef CUTCELLPROPERTIES_HPP_
#define CUTCELLPROPERTIES_HPP_

#include<iostream>
#include<string>

#include "TypeDefinition.hpp"
#include "geomPoint3D.hpp"
#include "geomPoint2D.hpp"
#include "delaunay.hpp"
#include "geomTriangle.hpp"
#include "geomTetra.hpp"
#include "geomHull.hpp"
#include "geomHull2D.hpp"
#include "geomCPgrid.hpp"
#include "geomFault.hpp"
#include "interGridIntersections.hpp"
#include "fracture.hpp"
#include "gmm/gmm.h"

extern "C" {
#ifdef _MSC_VER
# include <libqhull/qhull_a.h>
#else
# include <libqhull/libqhull.h>
  //# include <qhull/mem.h>
# include <libqhull/qset.h>
# include <libqhull/geom.h>
# include <libqhull/merge.h>
# include <libqhull/poly.h>
# include <libqhull/io.h>
# include <libqhull/stat.h>
#endif
}
namespace Geometry
{

  /*!
    @class Hull

    @author Anna Scotti

      */
  class CProp
  {
  public:

    typedef std::vector<Real> vettReal;
    typedef std::vector<Point3D> vettPoints;
    //! @name Constructor & Destructor
    //@{

    //! Empty constructor
    CProp();

    //! Constructor

    CProp (const Intersect::GridIntersections&, CPgrid*, Fracture*);

    //! Destructor
    virtual ~CProp();

    void setProperties();

    void setProperties_old();

    std::vector<Point3D> getIntPoints (Intersect::GridIntersections_Const_Iterator_Type&);

    std::vector<bool> getIsIntPointReal (Intersect::GridIntersections_Const_Iterator_Type&);

    void buildIntSegments (Intersect::GridIntersections_Const_Iterator_Type& it);

    void getCellPoints (Intersect::GridIntersections_Const_Iterator_Type&, std::vector<Point3D>&, int , Point3D const&)  const;

    Real setIntArea (Hull const& , gmm::size_type ) const;
    //! Computes the area of a grid of triangles
    /*!
    @param points Vector of points
    @param elements Vector of integere triplets identifing the triangular elements.
    @return The sum of the area of all triangles.
    */
    Real setIntArea (std::vector<IdTriplet> const& elements, std::vector<Point3D> const& points) const;

    Point3D setCG (std::vector<Point3D> const&) const;

    Point3D setCG (Hull const& , gmm::size_type ) const;

    Real setIntd (Hull const& , Point3D const&) const;

    //! Find  points for the computation of area and baricenter of the intersection.
    /*!
      It implements an algorithm that generates a set of points on the boundary of the fracture
      with the purpose of calculating the area of the intersecting region and its center.
      @param puntiarea The intersection points (real and virtual) found by the Newton algorithm.
      @param puntiIsReal True if the intersection point is a real one.
      @param it Iterator to the intersection
     */
    std::vector<Point3D> addPoints4area (std::vector<Point3D> const& puntiarea, std::vector<bool> const& puntiIsReal, Intersect::GridIntersections_Const_Iterator_Type& it) const;

    //! Find  points for the computation of area and baricenter of the intersection.
    /*!
      It implements an algorithm that generates a set of points on the boundary of the fracture
      with the purpose of calculating the area of the intersecting region and its center.
      @param puntiarea The intersection points (real and virtual) found by the Newton algorithm.
      @param puntiIsReal True if the intersection point is a real one.
      @param puntiuv Parametric coordinates of the point.
      @param it Iterator to the intersection.
     */
    std::vector<Point3D> addPoints4area (std::vector<Point3D> const& puntiarea, std::vector<bool> const& puntiIsReal, Intersect::GridIntersections_Const_Iterator_Type& it, std::vector<Point2D>& puntiuv) const;

    Real setIntdist_linea (std::vector<IdTriplet> const&, std::vector<Point3D> const& , Fracture::IntFrac const&) const;

    Real setIntdist_linea (Hull const&, gmm::size_type, Fracture::IntFrac const&) const;

    Real setIntdist_linea (std::vector<Point3D> const&, Point3D const&, Fracture::IntFrac const& , Intersect::GridIntersections_Const_Iterator_Type&, bool) const ;

    inline CPgrid* getgridpointer()
    {
      return M_gridpointer;
    }

    inline std::vector<Real>& getVolumes()
    {
      return M_vol;
    }

    inline std::vector<Real>& getAreas()
    {
      return M_aree;
    }

    inline std::vector<Real>& getDmedio()
    {
      return M_dmedio;
    }

    inline std::vector<Point3D>& getCG()
    {
      return M_CG;
    }

    inline std::vector<Segment>& getSx (gmm::size_type i)
    {
      return (i == 1) ? M_S1x : M_S2x;
    }
    inline std::vector<Segment>& getSy (gmm::size_type i)
    {
      return (i == 1) ? M_S1y : M_S2y;
    }
    inline std::vector<Segment>& getSz (gmm::size_type i)
    {
      return (i == 1) ? M_S1z : M_S2z;
    }

    inline gmm::size_type getNe() const
    {
      return M_Ne;
    }

    inline std::vector<Real> getMdmedioInt (gmm::size_type i) const
    {
      return M_dmedioint[i];
    }
    inline vettPoints getPoints (gmm::size_type i) const
    {
      return M_puntiAree[i];
    }


    inline std::vector<gmm::size_type> getI() const
    {
      return M_i;
    }
    inline std::vector<gmm::size_type> getJ() const
    {
      return M_j;
    }
    inline std::vector<gmm::size_type> getK() const
    {
      return M_k;
    }
    inline UInt getNx() const
    {
      return M_gridpointer->Nx();
    }
    inline UInt getNy() const
    {
      return M_gridpointer->Ny();
    }
    inline UInt getNz() const
    {
      return M_gridpointer->Nz();
    }
  private:
    std::vector<Real> M_vol, M_aree,  M_dmedio;
    std::vector<vettReal> M_dmedioint;
    std::vector<vettPoints> M_puntiAree;
    std::vector<Point3D> M_CG;
    Intersect::GridIntersections_Const_Iterator_Type M_iteratorcellsbegin, M_iteratorcellsend;
    CPgrid* M_gridpointer;
    Fracture* M_faultpointer;
    gmm::size_type M_Ne;
    std::vector<gmm::size_type> M_i;
    std::vector<gmm::size_type> M_j;
    std::vector<gmm::size_type> M_k;
    std::vector<Segment> M_S1x;
    std::vector<Segment> M_S1y;
    std::vector<Segment> M_S1z;
    std::vector<Segment> M_S2x;
    std::vector<Segment> M_S2y;
    std::vector<Segment> M_S2z;

  };



} // namespace Geometry

#endif /* GEOMTRIANGLE_HPP_ */
