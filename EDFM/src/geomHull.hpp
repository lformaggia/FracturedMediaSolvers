/*!
 *  @file geomTriangle.hpp
 *  @brief Class for Triangle in 3D space.
 *
 *  @author Luca Turconi <lturconi@gmail.com>
 *  @date 22-09-2012
 *
 */

#ifndef GEOMHULL_HPP_
#define GEOMHULL_HPP_

#include<iostream>
#include<string>

#include "TypeDefinition.hpp"
#include "geomPoint3D.hpp"
#include "geomTriangle.hpp"
#include "geomTetra.hpp"
#include "gmm/gmm.h"

extern "C" {
#ifdef _MSC_VER
# include <libqhull/qhull_a.h>
#else
# include <libqhull/libqhull.h>
# include <libqhull/mem.h>
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

    @brief It implements the interface to the QHull library.

    @author Anna Scotti

      */
  class Hull
  {
  public:

    //! @name Constructor & Destructor
    //@{

    //! Empty constructor
    Hull();

    //! Constructor, getting the simplexes matrix

    Hull (const gmm::dense_matrix<gmm::size_type>, const std::vector<coordT>);

    Hull (std::vector<coordT>);

    Hull (std::vector<Point3D>);

    //! Destructor
    virtual ~Hull();

    //@}

    //! @name Get Methods
    //@{

    //! Get vertex A
    /*!
     * @return The vertex A
     */
    inline gmm::size_type getNtetra() const
    {
      return M_Ntetra;
    }
    inline Point3D getPoint (gmm::size_type& i) const
    {
      return M_points[i];
    }
    inline gmm::dense_matrix<gmm::size_type> getSimplexes() const
    {
      return M_simplexes;
    }

    bool call_qhull (gmm::size_type& , gmm::dense_matrix<gmm::size_type>&, std::vector<coordT>&);

    inline std::vector<Tetra> getTetra() const
    {
      return M_tetra;
    }

    inline Tetra getTetra (gmm::size_type i) const
    {
      return M_tetra[i];
    }

    Real getVolume() const;

    std::vector<gmm::size_type> getPointsSimplex (gmm::size_type ) const;

  private:
    gmm::dense_matrix<gmm::size_type> M_simplexes;
    std::vector<Point3D> M_points;
    std::vector<Tetra> M_tetra;
    gmm::size_type M_Ntetra;
  };



} // namespace Geometry

#endif /* GEOMTRIANGLE_HPP_ */
