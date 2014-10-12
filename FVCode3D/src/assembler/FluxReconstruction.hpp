/*!
 *  @file FluxReconstruction.hpp
 *  @brief This class computes the flux/velocity from the Darcy problem.
 */

#ifndef __FLUX_HPP__
#define __FLUX_HPP__

#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
#include "boundaryCondition/BC.hpp"
#include "core/TypeDefinition.hpp"
#include "mesh/Rigid_Mesh.hpp"

class Rigid_mesh;
class PropertiesMap;

namespace Darcy
{

//! Class for reconstruction of the flux
/*!
    @class FluxReconstruction
    This class computes the flux/velocity from the Darcy problem.
*/
class FluxReconstruction: public StiffMatrix
{
    
    //! Typedef for std::pair<UInt,UInt>
    /*!
        @typedef Fracture_Juncture
        This type definition permits to treat std::pair<UInt,UInt> as a Fracture_Juncture.
    */
    typedef std::pair<UInt,UInt> Fracture_Juncture;
    //! Typedef for Geometry::Rigid_Mesh::Facet_ID
    /*!
        @typedef Facet_ID
        This type definition permits to treat Geometry::Rigid_Mesh::Facet_ID as a Facet_ID.
    */
    typedef Geometry::Rigid_Mesh::Facet_ID Facet_ID;
    
public:

    //! @name Constructor & Destructor
    //@{
    //! Constructor
    /*!
     * @param pressure vector of the pressure on cells
     */ 
    FluxReconstruction(const Geometry::Rigid_Mesh & rigid_mesh, const BoundaryConditions & bc, const std::vector<Real> & pressure):
        M_mesh(rigid_mesh), M_bc(bc), M_properties(rigid_mesh.getPropertiesMap()), M_pressure(pressure) {};
        
    //! No Copy-Constructor
    FluxReconstruction(const FluxReconstruction &) = delete;
    
    //! No Empty-Constructor
    FluxReconstruction() = delete;

    //! Destructor
    virtual ~FluxReconstruction() {};
    //@}

    //! @name Get Methods
    //@{    
    const std::vector<Real> & getFlux() const
        { return M_flux; };
    
    std::vector<Real> & getFlux()
        { return M_flux; };
        
    const std::vector<Real> & getVelocity() const
        { return M_velocity; };
        
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
    
    //! A reference to a Geometry::Rigid_Mesh
    const Geometry::Rigid_Mesh & M_mesh;
    //! A reference to a Geometry::PropertiesMap
    const Geometry::PropertiesMap & M_properties;
    //! The container of the BCs
    const BoundaryConditions & M_bc;

};

} // namespace Darcy

#endif