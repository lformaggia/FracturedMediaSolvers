/*!
 * @file exportVTK.cpp
 * @brief Classes for saving files in the VTK format (definitions).
 */
#include "export/exportVTK.hpp"

void ExporterVTK::exportMesh(const Mesh3D & mesh, const std::string filename)
{
    // TODO
}

void ExporterVTK::exportFractures(const Mesh3D & mesh, const std::string filename)
{
    // TODO
}

void ExporterVTK::exportMeshWithFractures(const Mesh3D & mesh, const std::string filename)
{
    // TODO
}

void ExporterVTK::exportFractureJunctures(const Rigid_Mesh & mesh, const std::string filename)
{
    // TODO
}

void ExporterVTK::exportSolutionOnFractures(const Rigid_Mesh & mesh, const std::string filename, const Eigen::VectorXd & sol)
{
    // TODO
}

void ExporterVTK::exportWithProperties(const Mesh3D & mesh, const PropertiesMap & properties, const std::string filename, const Flag16bit propertiesType, const std::vector<Real> * property)
{
    // TODO
}

void ExporterVTK::exportWithProperties(const Mesh3D & mesh, const PropertiesMap & properties, const std::string filename)
{
    // TODO
}

void ExporterVTK::exportWithProperties(const Rigid_Mesh & mesh, const std::string filename, const Flag16bit propertiesType, const std::vector<Real> * property)
{
    // TODO
}
