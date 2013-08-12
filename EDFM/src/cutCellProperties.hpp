/*!
 *  @file geomTriangle.hpp
 *  @brief Class for Triangle in 3D space.
 *
 *  @author Luca Turconi <lturconi@gmail.com>
 *  @date 22-09-2012
 *
 */

#ifndef CUTCELLPROPERTIES_HPP_
#define CUTCELLPROPERTIES_HPP_

#include<iostream>
#include<string>

#include "TypeDefinition.hpp"
#include "geomPoint3D.hpp"
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

    std::vector<Point3D> getIntPoints (Intersect::GridIntersections_Const_Iterator_Type&);

    std::vector<bool> getIsIntPointReal (Intersect::GridIntersections_Const_Iterator_Type&);

    void buildIntSegments (Intersect::GridIntersections_Const_Iterator_Type& it);

    void getCellPoints (Intersect::GridIntersections_Const_Iterator_Type&, std::vector<Point3D>&, int , Point3D const &)  const;

    Real setIntArea (Hull const & , gmm::size_type ) const;

    Point3D setCG (std::vector<Point3D> const &)const;

    Point3D setCG (Hull const & , gmm::size_type )const;

    Real setIntd (Hull& , Point3D);

    std::vector<Point3D> addPoints4area (std::vector<Point3D> const&, std::vector<bool> const &, Intersect::GridIntersections_Const_Iterator_Type&) const;

    Real setIntdist_linea (Hull&, gmm::size_type, Fracture::IntFrac);

    Real setIntdist_linea (std::vector<Point3D>, Point3D, Fracture::IntFrac, Intersect::GridIntersections_Const_Iterator_Type&, bool);

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

    inline std::vector<Real> getMdmedioInt (gmm::size_type i)
    {
      return M_dmedioint[i];
    }
    inline vettPoints getPoints (gmm::size_type i)
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
