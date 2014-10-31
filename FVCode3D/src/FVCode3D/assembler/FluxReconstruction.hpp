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
     * @param pressure vector of the pressure on cells
     */
    FluxReconstruction(const Rigid_Mesh & rigid_mesh, const BoundaryConditions & bc, const std::vector<Real> & pressure):
        StiffMatrix(rigid_mesh, bc), M_pressure(pressure)
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
    const std::vector<Real> & getFlux() const
    { return M_flux; };

    //! Get the flux
    /*!
     * @return the flux vector
     */
    std::vector<Real> & getFlux()
    { return M_flux; };

    //! Get the velocity (const)
    /*!
     * @return the velocity vector
     */
    const std::vector<Real> & getVelocity() const
    { return M_velocity; };

    //! Get the velocity (const)
    /*!
     * @return the velocity vector
     */
    std::vector<Real> & getVelocity()
    { return M_velocity; };
    //@}

    //! @name Methods
    //@{
    //! Execute the reconstruction of flux/velocity
    void reconstruct();
    //@}

private:

    //! Vector that contains the flux on the facets
    std::vector<Real> M_flux;
    //! Vector that contains the velocity on each cell: [ Vx1 ... VxN  Vy1 ... VyN  Vz1 ... VzN ]
    std::vector<Real> M_velocity;

    //! Vector that contains the pressure
    const std::vector<Real> & M_pressure;
}; // namespace FluxReconstruction

} // namespace FVCode3D

#endif // __FLUX_HPP__
