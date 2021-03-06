/*!
 * @file export.hpp
 * @brief Classes for saving files.
 */

#ifndef EXPORT_HPP_
#define EXPORT_HPP_

#include <FVCode3D/core/TypeDefinition.hpp>
#include <FVCode3D/mesh/RigidMesh.hpp>

namespace FVCode3D
{

class Mesh3D;
class Rigid_Mesh;
class PropertiesMap;

//! Flag for exporter
/*!
 * @enum PropertiesToExport
 * This enum permits to select which properties to export.
 */
enum PropertiesToExport
{
    ZoneCode        = 0x001,
    IsFrac          = 0x002,
    Aperture        = 0x004,
    Porosity        = 0x008,
    Permeability    = 0x010,
    BorderID        = 0x020,
    ElementID       = 0x040,
    Other           = 0x080,
};

//! Class to export mesh, fractures, solution and properties.
/*!
 * @class Exporter
 * This is a base abstact class that allows to export mesh, fractures, solution, and properties.
 * Each derived class permits to export in a specific file format.
 */
class Exporter
{
public:

    typedef Mesh3D::Cell3D Cell3D;
    typedef Mesh3D::Facet3D Facet3D;

    typedef Rigid_Mesh::Cell Cell;
    typedef Rigid_Mesh::Facet Facet;
    typedef Rigid_Mesh::Edge Edge;

    typedef Rigid_Mesh::Fracture_Facet Fracture_Facet;
    typedef Rigid_Mesh::Fracture_Tip Fracture_Tip;
    typedef Rigid_Mesh::Fracture_Juncture Fracture_Juncture;

    typedef Rigid_Mesh::Regular_Edge Regular_Edge;
    typedef Rigid_Mesh::Juncture_Edge Juncture_Edge;
    typedef Rigid_Mesh::Internal_Tip_Edge Internal_Tip_Edge;
    typedef Rigid_Mesh::Border_Tip_Edge Border_Tip_Edge;
    typedef Rigid_Mesh::Pure_Border_Edge Pure_Border_Edge;

    //! Constructor
    Exporter() = default;

    //! Export the mesh (only cells)
    /*!
     * @param mesh reference of a Mesh3D
     * @param filename name of the file
     */
    virtual void exportMesh(const Mesh3D & mesh, const std::string filename) const= 0;

    //! Export a tetrahedral mesh (only cells)
    /*!
     * @param mesh reference of a Mesh3D
     * @param filename name of the file
     */
    virtual void exportTetrahedralMesh(const Mesh3D & mesh, const std::string filename) const = 0;

    //! Export the fracture facets
    /*!
     * @param mesh reference of a Mesh3D
     * @param filename name of the file
     */
    virtual void exportFractures(const Mesh3D & mesh, const std::string filename) const = 0;

    //! Export the mesh, cells and fracture facets, in a single file
    /*!
     * @param mesh reference of a Mesh3D
     * @param filename name of the file
     */
    virtual void exportMeshWithFractures(const Mesh3D & mesh, const std::string filename) const = 0;

    //! Export the wireframe
    /*!
     * @param mesh reference of a Mesh3D
     * @param filename name of the file
     */
    virtual void exportWireframe(const Mesh3D & mesh, const std::string filename) const = 0;

    //! Export the edges
    /*!
     * @param mesh reference of a Rigid_Mesh
     * @param filename name of the file
     */
    virtual void exportEdges(const Rigid_Mesh & mesh, const std::string filename) const = 0;

    //! Export the faces
    /*!
     * @param mesh reference of a Rigid_Mesh
     * @param filename name of the file
     */
    virtual void exportFacets(const Rigid_Mesh & mesh, const std::string filename) const = 0;

    //! Export the fracture junctures
    /*!
     * @param mesh reference of a Rigid_Mesh
     * @param filename name of the file
     */
    virtual void exportFractureJunctures(const Rigid_Mesh & mesh, const std::string filename) const = 0;

    //! Export the fracture tips
    /*!
     * @param mesh reference of a Rigid_Mesh
     * @param filename name of the file
     */
    virtual void exportFractureTips(const Rigid_Mesh & mesh, const std::string filename) const = 0;

    //! Export the solution on cells and fracture facets in a single file
    /*!
     * @tparam VectorType vector type (e.g., std::vector, Eigen::Vector ...)
     * @param mesh reference of a Rigid_Mesh
     * @param filename name of the file
     * @param sol Eigen vector that contain the solution (cells + fracture facets)
     * @param fieldName name of the field that appears in the file
     */
    template <typename VectorType>
    void exportSolution(const Rigid_Mesh & mesh, const std::string filename, const VectorType & sol, const std::string & fieldName = "Pressure") const;

    //! Export the solution on fracture facets
    /*!
     * @tparam VectorType vector type (e.g., std::vector, Eigen::Vector ...)
     * @param mesh reference of a Rigid_Mesh
     * @param filename name of the file
     * @param sol Eigen vector that contain the solution (cells + fracture facets)
     * @param fieldName name of the field that appears in the file
     */
    template <typename VectorType>
    void exportSolutionOnFractures(const Rigid_Mesh & mesh, const std::string filename, const VectorType & sol, const std::string & fieldName = "Pressure") const;

    //! Export the solution on cells and fracture facets in a single file
    /*!
     * @tparam VectorType vector type (e.g., std::vector, Eigen::Vector ...)
     * @param mesh reference of a Rigid_Mesh
     * @param filename name of the file
     * @param sol Eigen vector that contain the solution (cells + fracture facets)
     * @param fieldName name of the field that appears in the file
     */
    template <typename VectorType>
    void exportFlux(const Rigid_Mesh & mesh, const std::string filename, const VectorType & sol, const std::string & fieldName = "Flux") const;

    //! Export the solution on fracture facets
    /*!
     * @tparam VectorType vector type (e.g., std::vector, Eigen::Vector ...)
     * @param mesh reference of a Geometry::Rigid_Mesh
     * @param filename name of the file
     * @param sol Eigen vector that contain the solution (cells + fracture facets)
     * @param fieldName name of the field that appears in the file
     */
    template <typename VectorType>
    void exportFluxOnFractures(const Rigid_Mesh & mesh, const std::string filename, const Eigen::VectorXd & sol, const std::string & fieldName = "Flux") const;

    //! Export a solution on facets
    /*!
     * @param mesh reference of a Rigid_Mesh
     * @param filename name of the file
     * @param sol a vector to be exported
     * @param fieldName name of the field that appears in the file
     */
    template <typename VectorType>
    void exportSolutionOnFacets(const Rigid_Mesh & mesh, const std::string filename, const VectorType & sol, const std::string & fieldName = "Transmissibility") const;

    //! Export the a specific property on cells and fracture facets
    /*!
     * @param mesh reference of a Mesh3D
     * @param properties reference to a PropertiesMap
     * @param filename name of the file
     * @param propertiesType flag used to select which properties to export
     * @param property pointer to a generic property
     */
    virtual void exportWithProperties(const Mesh3D & mesh, const PropertiesMap & properties, const std::string filename,
        const Flag16bit propertiesType, const std::vector<Real> * property = nullptr ) const= 0;

    //! Export all properties defined on cells and fracture facets
    /*!
     * @param mesh reference of a Mesh3D
     * @param properties reference to a PropertiesMap
     * @param filename name of the file
     */
    virtual void exportWithProperties(const Mesh3D & mesh, const PropertiesMap & properties, const std::string filename) const = 0;

    //! Export the solution on fracture facets
    /*!
     * @param mesh reference of a Rigid_Mesh
     * @param filename name of the file
     * @param sol Eigen vector that contain the solution (cells + fracture facets)
     * @param property pointer to a generic property
     */
    virtual void exportWithProperties(const Rigid_Mesh & mesh, const std::string filename,
        const Flag16bit propertiesType, const std::vector<Real> * property = nullptr ) const= 0;

    //! Destructor
    virtual ~Exporter() = default;

};

} // namespace FVCode3D
#endif /* EXPORT_HPP_ */
