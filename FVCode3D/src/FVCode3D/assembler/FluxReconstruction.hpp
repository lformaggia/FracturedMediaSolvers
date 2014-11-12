/*!
 * @file FluxReconstruction.hpp
 * @brief This class computes the flux/velocity from the Darcy problem.
 */

#ifndef __FLUX_HPP__
#define __FLUX_HPP__

#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
#include <FVCode3D/assembler/Stiffness.hpp>

namespace FVCode3D
{

class Rigid_mesh;
class PropertiesMap;

//! Class for reconstruction of the flux
/*!
 * @class FluxReconstruction
 * This class computes the flux/velocity from the Darcy problem.
 */
class FluxReconstruction: public StiffMatrix
{

    //! Typedef for std::pair<UInt,UInt>
    /*!
     * @typedef Fracture_Juncture
     * This type definition permits to treat std::pair<UInt,UInt> as a Fracture_Juncture.
     */
    typedef std::pair<UInt,UInt> Fracture_Juncture;

    //! Typedef for Geometry::Rigid_Mesh::Facet_ID
    /*!
     * @typedef Facet_ID
     * This type definition permits to treat Geometry::Rigid_Mesh::Facet_ID as a Facet_ID.
     */
    typedef Rigid_Mesh::Facet_ID Facet_ID;

public:

    //! @name Constructor & Destructor
    //@{
    //! Constructor
    /*!
     * @param rigid_mesh A Rigid_Mesh used to build the matrix
     * @param BC Boundary conditions given in the container BoundaryConditions
     */
    FluxReconstruction(const Rigid_Mesh & rigid_mesh, const BoundaryConditions & bc):
        StiffMatrix(rigid_mesh, bc)
    {};

    //! No Copy-Constructor
    FluxReconstruction(const FluxReconstruction &) = delete;

    //! No Empty-Constructor
    FluxReconstruction() = delete;

    //! Destructor
    virtual ~FluxReconstruction() = default;
    //@}

    //! @name Get Methods
    //@{
    //! Get the flux (const)
    /*!
     * @return the flux vector
     */
    const Vector & getFlux() const
    { return M_flux; };

    //! Get the flux
    /*!
     * @return the flux vector
     */
    Vector & getFlux()
    { return M_flux; };
    //@}

    //! @name Methods
    //@{
    //! Execute the reconstruction of flux
    void reconstruct(const Vector & pressure);

    //! Execute the reconstruction of flux
    void localReconstruct(const UInt facetId, const Vector & pressure);
    //@}

private:
    //! Vector that contains the flux on the facets
    Vector M_flux;
}; // FluxReconstruction

} // namespace FVCode3D

#endif // __FLUX_HPP__
