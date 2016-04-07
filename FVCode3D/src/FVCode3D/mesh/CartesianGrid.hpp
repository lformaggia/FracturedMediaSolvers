/*!
 * @file CartesianGrid.hpp
 * @brief Class that generate a hexahedral structured (Cartesian) grid.
 */

#ifndef CARTESIANGRID_HPP_
#define CARTESIANGRID_HPP_

#include <functional>

#include <FVCode3D/core/TypeDefinition.hpp>
#include <FVCode3D/core/Data.hpp>

namespace FVCode3D
{

class Mesh3D;
class PropertiesMap;

//! Class used to generate hexahedaral structured meshes.
/*!
 * @class CartesianGrid
 * This class that allows to generate Cartesian grids.
 */
class CartesianGrid
{
public:

    //! Constructor
    /*!
     * @param mesh reference to a Mesh3D
     * @param properties reference to a PropertiesMap
     * @param data reference to a Data class
     */
    CartesianGrid(Mesh3D & mesh, PropertiesMap & properties, const DataPtr_Type & data):
        M_mesh(mesh), M_properties(properties), M_data(data) {}

    //! Generate a Cartesian mesh with a prescribed surface
    /*!
     * @param fracturesOn true to load fractures
     */
    void generate( bool fracturesOn,
                   const std::function<Real(const Point3D&)> & _surface );

    //! Generate a Cartesian mesh
    /*!
     * @param fracturesOn true to load fractures
     */
    virtual void generate( bool fracturesOn )
    { generate( fracturesOn,
                [&](const Point3D&){ return M_data->getLz(); } ); }

    //! Generate the BC ids
    /*!
     * Add the BCs ids to the boundary facets.
     * The id is set by considering the maximum component and the sign of the normal of a boundary facet.
     * @param theta rotation angle along z-axis. It is used only to compute the BC ids. Default = 0
     * @pre generate()
     */
    virtual void extractBC(const Real theta = 0.);

    //! Add the fractures network
    /*!
     * Add the fracture network from a list of facet ids
     * @param facetIdToZone map that associates to a facet id the zone code
     * @pre generate()
     */
    virtual void addFractures(const std::map<UInt,UInt> & facetIdToZone);

    //! Add the lacking data
    /*!
     * Add the BCs and the fracture network
     * @param facetIdToZone map that associates to a facet id the zone code
     * @param theta rotation angle along z-axis. It is used only to compute the BC ids. Default = 0
     * @pre generate()
     */
    virtual void addBCAndFractures(const std::map<UInt,UInt> & facetIdToZone, const Real theta = 0.);

    //! Get mesh
    /*!
     * @return reference to the mesh
     */
    const Mesh3D & getMesh() const
        { return M_mesh; }

    //! Get properties
    /*!
     * @return reference to the properties
     */
    const PropertiesMap & getProperties() const
        { return M_properties; }

    //! Destructor
    virtual ~CartesianGrid() = default;

protected:

    //! Reference to a Mesh3D
    Mesh3D & M_mesh;
    //! Reference to a PropertiesMap
    PropertiesMap & M_properties;
    //! Pointer to Data
    const DataPtr_Type & M_data;

private:

    //! No default constructor
    CartesianGrid() = delete;

    //! No copy-constructor
    CartesianGrid(const CartesianGrid &) = delete;

    //! No assignment operator
    CartesianGrid & operator=(const CartesianGrid &) = delete;
};

} // FVCode3D
#endif /* CARTESIANGRID_HPP_ */
