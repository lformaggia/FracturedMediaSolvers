/*
 * @file ExportCP.hpp
 * @brief Classes for saving files in the CP format (GRDECL).
 */

#ifndef EXPORTCP_HPP_
#define EXPORTCP_HPP_

#include <FVCode3D/export/Export.hpp>

namespace FVCode3D
{

/*!
 * @class ExporterCP
 * This class allows to export mesh and fractures.
 */
class ExporterCP : public Exporter
{
public:

    //! Constructor
    ExporterCP() = default;

    //! Export the mesh (only cells)
    /*!
     * @param mesh reference of a Geometry::Mesh3D
     * @param filename name of the file
     */
    virtual void exportMesh(const Mesh3D & mesh, const std::string filename) throw();

    //! Destructor
    virtual ~ExporterCP() = default;

protected:

    //! Export the fracture facets
    /*!
     * @param mesh reference of a Geometry::Mesh3D
     * @param filename name of the file
     */
    virtual void exportFractures(const Mesh3D & /*mesh*/, const std::string /*filename*/) throw(){}

    //! Export a tetrahedral mesh (only cells)
    /*!
     * @param mesh reference of a Geometry::Mesh3D
     * @param filename name of the file
     */
    virtual void exportTetrahedralMesh(const Mesh3D & /*mesh*/, const std::string /*filename*/) throw(){}

    //! Export the mesh, cells and fracture facets, in a single file
    /*!
     * @param mesh reference of a Geometry::Mesh3D
     * @param filename name of the file
     */
    virtual void exportMeshWithFractures(const Mesh3D & /*mesh*/, const std::string /*filename*/) throw(){}

    //! Export the wireframe
    /*!
     * @param mesh reference of a Geometry::Mesh3D
     * @param filename name of the file
     */
    virtual void exportWireframe(const Mesh3D & /*mesh*/, const std::string /*filename*/) throw(){}

    //! Export the edges
    /*!
     * @param mesh reference of a Rigid_Mesh
     * @param filename name of the file
     */
    virtual void exportEdges(const Rigid_Mesh & /*mesh*/, const std::string /*filename*/) throw(){}

    //! Export the faces
    /*!
     * @param mesh reference of a Rigid_Mesh
     * @param filename name of the file
     */
    virtual void exportFacets(const Rigid_Mesh & /*mesh*/, const std::string /*filename*/) throw(){}

    //! Export the fracture junctures
    /*!
     * @param mesh reference of a Geometry::Rigid_Mesh
     * @param filename name of the file
     */
    virtual void exportFractureJunctures(const Rigid_Mesh & /*mesh*/, const std::string /*filename*/) throw(){}

    //! Export the fracture tips
    /*!
     * @param mesh reference of a Geometry::Rigid_Mesh
     * @param filename name of the file
     */
    virtual void exportFractureTips(const Rigid_Mesh & /*mesh*/, const std::string /*filename*/) throw(){}

    //! Export the solution on cells and fracture facets in a single file
    /*!
     * @param mesh reference of a Geometry::Rigid_Mesh
     * @param filename name of the file
     * @param sol Eigen vector that contain the solution (cells + fracture facets)
     */
    template <typename VectorType>
    void exportSolution(const Rigid_Mesh & /*mesh*/, const std::string /*filename*/, const VectorType & /*sol*/, const std::string & /*fieldName = "Pressure"*/) throw(){}

    //! Export the solution on fracture facets
    /*!
     * @param mesh reference of a Geometry::Rigid_Mesh
     * @param filename name of the file
     * @param sol Eigen vector that contain the solution (cells + fracture facets)
     */
    template <typename VectorType>
    void exportSolutionOnFractures(const Rigid_Mesh & /*mesh*/, const std::string /*filename*/, const VectorType & /*sol*/, const std::string & /*fieldName = "Pressure"*/) throw(){}

    //! Export the a specific property on cells and fracture facets
    /*!
     * @param mesh reference of a Geometry::Mesh3D
     * @param properties reference to a Geometry::PropertiesMap
     * @param filename name of the file
     * @param propertiesType flag used to select which properties to export
     * @param property pointer to a generic property
     */
    virtual void exportWithProperties(const Mesh3D & /*mesh*/, const PropertiesMap & /*properties*/, const std::string /*filename*/, const Flag16bit /*propertiesType*/, const std::vector<Real> * /*property = nullptr*/) throw(){}

    //! Export all properties defined on cells and fracture facets
    /*!
     * @param mesh reference of a Geometry::Mesh3D
     * @param properties reference to a Geometry::PropertiesMap
     * @param filename name of the file
     */
    virtual void exportWithProperties(const Mesh3D & /*mesh*/, const PropertiesMap & /*properties*/, const std::string /*filename*/) throw(){}

    //! Export the solution on fracture facets
    /*!
     * @param mesh reference of a Geometry::Rigid_Mesh
     * @param filename name of the file
     * @param sol Eigen vector that contain the solution (cells + fracture facets)
     */
    virtual void exportWithProperties(const Rigid_Mesh & /*mesh*/, const std::string /*filename*/, const Flag16bit /*propertiesType*/, const std::vector<Real> * /*property = nullptr*/) throw(){}
};

} // namespace FVCode3D

#endif /* EXPORTCP_HPP_ */
