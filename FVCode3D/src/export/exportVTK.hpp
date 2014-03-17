/*!
 *  @file exportVTK.hpp
 *  @brief Classes for saving files in the VTK format.
 */

#ifndef EXPORTVTK_HPP_
#define EXPORTVTK_HPP_

#include "export/export.hpp"

//! Class to export mesh, fractures, solution and properties in VTK format
/*!
 * @class ExporterVTK
 * This class allows to export mesh, fractures, solution, and properties in VTK format.
 * @warning Not implemented.
 */
class ExporterVTK : public Exporter
{
public:

    ExporterVTK(){}

    virtual void exportMesh(const Mesh3D & mesh, const std::string filename);

    virtual void exportFractures(const Mesh3D & mesh, const std::string filename);

    virtual void exportMeshWithFractures(const Mesh3D & mesh, const std::string filename);

    virtual void exportFractureJunctures(const Rigid_Mesh & mesh, const std::string filename);

    template <typename VectorType>
    void exportSolution(const Rigid_Mesh & mesh, const std::string filename, const VectorType & sol);

    virtual void exportSolutionOnFractures(const Rigid_Mesh & mesh, const std::string filename, const Eigen::VectorXd & sol);

    virtual void exportWithProperties(const Mesh3D & mesh, const PropertiesMap & properties, const std::string filename, const Flag16bit propertiesType, const std::vector<Real> * property = NULL);

    virtual void exportWithProperties(const Mesh3D & mesh, const PropertiesMap & properties, const std::string filename);

    virtual void exportWithProperties(const Rigid_Mesh & mesh, const std::string filename, const Flag16bit propertiesType, const std::vector<Real> * property = NULL );

    virtual ~ExporterVTK() {};
};

template <typename VectorType>
void ExporterVTK::exportSolution(const Rigid_Mesh & mesh, const std::string filename, const VectorType & sol)
{
    // TODO
}

#endif /* EXPORTVTK_HPP_ */
