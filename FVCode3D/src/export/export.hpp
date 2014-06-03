/*!
 *  @file export.hpp
 *  @brief Classes for saving files.
 */

#ifndef EXPORT_HPP_
#define EXPORT_HPP_

#include "core/TypeDefinition.hpp"
#include "mesh/Rigid_Mesh.hpp"

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

    typedef Geometry::Mesh3D Mesh3D;
    typedef Geometry::Mesh3D::Cell3D Cell3D;
    typedef Geometry::Mesh3D::Facet3D Facet3D;
    typedef Geometry::PropertiesMap PropertiesMap;

    typedef Geometry::Rigid_Mesh Rigid_Mesh;
    typedef Geometry::Rigid_Mesh::Cell Cell;
    typedef Geometry::Rigid_Mesh::Facet Facet;
    typedef Geometry::Rigid_Mesh::Fracture_Facet Fracture_Facet;
    typedef Geometry::Rigid_Mesh::Fracture_Tip Fracture_Tip;
    typedef Geometry::Rigid_Mesh::Fracture_Juncture Fracture_Juncture;

    //! Constructor
    Exporter(){}

    //! Export the mesh (only cells)
    /*!
     * @param mesh reference of a Geometry::Mesh3D
     * @param filename name of the file
     */
    virtual void exportMesh(const Mesh3D & mesh, const std::string filename) = 0;

    //! Export a tetrahedral mesh (only cells)
    /*!
     * @param mesh reference of a Geometry::Mesh3D
     * @param filename name of the file
     */
    virtual void exportTetrahedralMesh(const Mesh3D & mesh, const std::string filename) = 0;

    //! Export the fracture facets
    /*!
     * @param mesh reference of a Geometry::Mesh3D
     * @param filename name of the file
     */
    virtual void exportFractures(const Mesh3D & mesh, const std::string filename) = 0;

    //! Export the mesh, cells and fracture facets, in a single file
    /*!
     * @param mesh reference of a Geometry::Mesh3D
     * @param filename name of the file
     */
    virtual void exportMeshWithFractures(const Mesh3D & mesh, const std::string filename) = 0;

    //! Export the wireframe
    /*!
     * @param mesh reference of a Geometry::Mesh3D
     * @param filename name of the file
     */
    virtual void exportWireframe(const Mesh3D & mesh, const std::string filename) = 0;
    
    //! Export the fracture junctures
    /*!
     * @param mesh reference of a Geometry::Rigid_Mesh
     * @param filename name of the file
     */
    virtual void exportFractureJunctures(const Rigid_Mesh & mesh, const std::string filename) = 0;

    //! Export the fracture tips
    /*!
     * @param mesh reference of a Geometry::Rigid_Mesh
     * @param filename name of the file
     */
    virtual void exportFractureTips(const Rigid_Mesh & mesh, const std::string filename) = 0;

    //! Export the solution on cells and fracture facets in a single file
    /*!
     * @param mesh reference of a Geometry::Rigid_Mesh
     * @param filename name of the file
     * @param sol Eigen vector that contain the solution (cells + fracture facets)
     */
    template <typename VectorType>
    void exportSolution(const Rigid_Mesh & mesh, const std::string filename, const VectorType & sol);

    //! Export the solution on fracture facets
    /*!
     * @param mesh reference of a Geometry::Rigid_Mesh
     * @param filename name of the file
     * @param sol Eigen vector that contain the solution (cells + fracture facets)
     */
    virtual void exportSolutionOnFractures(const Rigid_Mesh & mesh, const std::string filename, const Eigen::VectorXd & sol) = 0;

    //! Export the solution on cells and fracture facets in a single file
    /*!
     * @param mesh reference of a Geometry::Rigid_Mesh
     * @param filename name of the file
     * @param sol Eigen vector that contain the solution (cells + fracture facets)
     */
    template <typename VectorType>
    void exportFlux(const Rigid_Mesh & mesh, const std::string filename, const VectorType & sol);

    //! Export the solution on fracture facets
    /*!
     * @param mesh reference of a Geometry::Rigid_Mesh
     * @param filename name of the file
     * @param sol Eigen vector that contain the solution (cells + fracture facets)
     */
    virtual void exportFluxOnFractures(const Rigid_Mesh & mesh, const std::string filename, const Eigen::VectorXd & sol) = 0;
    
    //! Export the a specific property on cells and fracture facets
    /*!
     * @param mesh reference of a Geometry::Mesh3D
     * @param properties reference to a Geometry::PropertiesMap
     * @param filename name of the file
     * @param propertiesType flag used to select which properties to export
     * @param property pointer to a generic property
     */
    virtual void exportWithProperties(const Mesh3D & mesh, const PropertiesMap & properties, const std::string filename, const Flag16bit propertiesType, const std::vector<Real> * property = NULL ) = 0;

    //! Export all properties defined on cells and fracture facets
    /*!
     * @param mesh reference of a Geometry::Mesh3D
     * @param properties reference to a Geometry::PropertiesMap
     * @param filename name of the file
     */
    virtual void exportWithProperties(const Mesh3D & mesh, const PropertiesMap & properties, const std::string filename) = 0;

    //! Export the solution on fracture facets
    /*!
     * @param mesh reference of a Geometry::Rigid_Mesh
     * @param filename name of the file
     * @param sol Eigen vector that contain the solution (cells + fracture facets)
     */
    virtual void exportWithProperties(const Rigid_Mesh & mesh, const std::string filename, const Flag16bit propertiesType, const std::vector<Real> * property = NULL ) = 0;

    //! Destructor
    virtual ~Exporter(){};

};

#endif /* EXPORT_HPP_ */
