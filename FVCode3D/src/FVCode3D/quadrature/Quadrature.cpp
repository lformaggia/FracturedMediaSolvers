/*!
 *  @file Quadrature.cpp
 *  @brief Class for integrating functions on a Rigid_Mesh (definitions).
 */

#include <FVCode3D/mesh/RigidMesh.hpp>
#include <FVCode3D/property/Properties.hpp>
#include <FVCode3D/quadrature/Quadrature.hpp>

namespace FVCode3D
{

Quadrature::Quadrature (const Rigid_Mesh & rigid_mesh):
     M_mesh (rigid_mesh), M_properties(rigid_mesh.getPropertiesMap()),
     M_size (rigid_mesh.getCellsVector().size()+ rigid_mesh.getFractureFacetsIdsVector().size()),
     M_quadrature(std::move(std::unique_ptr< QuadratureRule >(new CentroidQuadrature))),
     M_fractureQuadrature(std::move(std::unique_ptr< QuadratureRule >(new CentroidQuadrature)))
{}

Quadrature::Quadrature (const Rigid_Mesh & rigid_mesh, const QuadratureRule & quadrature):
    M_mesh (rigid_mesh), M_properties(rigid_mesh.getPropertiesMap()),
    M_size (rigid_mesh.getCellsVector().size()+ rigid_mesh.getFractureFacetsIdsVector().size()),
    M_quadrature(std::move(quadrature.clone())),
    M_fractureQuadrature(std::move(std::unique_ptr<QuadratureRule >(new CentroidQuadrature)))
{}

Quadrature::Quadrature (const Rigid_Mesh & rigid_mesh, const QuadratureRule & quadrature, const QuadratureRule & fracturequadrature):
    M_mesh (rigid_mesh), M_properties(rigid_mesh.getPropertiesMap()),
    M_size (rigid_mesh.getCellsVector().size()+ rigid_mesh.getFractureFacetsIdsVector().size()),
    M_quadrature(std::move(quadrature.clone())),
    M_fractureQuadrature(std::move(fracturequadrature.clone()))
{}

Real Quadrature::integrate(const Vector & integrand) const
{
    UInt IntSize = integrand.size();
    if(IntSize != M_size)
    {
        std::stringstream error;
        error << "Error: dimension of integrand function " << IntSize << " differs from dimension of mesh " << M_size << ".";
        throw std::runtime_error(error.str());
    }
    Real integral = 0;
    for(auto& cell_it : M_mesh.getCellsVector())
    {
        integral += cell_it.getVolume()*integrand(cell_it.getId());
    }

    for(auto& facet_it : M_mesh.getFractureFacetsIdsVector())
    {
//        Real _volume = M_properties.getProperties(facet_it.getZoneCode()).M_aperture * facet_it.getFacet().area();
        integral += facet_it.getFacet().area()*integrand(facet_it.getIdAsCell());
    }

    return integral;
}

Real Quadrature::integrateMatrix(const Vector & integrand) const
{
    UInt IntSize = integrand.size();
    if(IntSize != M_size)
    {
        std::stringstream error;
        error << "Error: dimension of integrand function " << IntSize << " differs from dimension of mesh " << M_size << ".";
        throw std::runtime_error(error.str());
    }
    Real integral = 0;
    for(auto& cell_it : M_mesh.getCellsVector())
    {
        integral += cell_it.getVolume()*integrand(cell_it.getId());
    }
    return integral;
}

Real Quadrature::integrateFractures(const Vector & integrand) const
{
    UInt IntSize = integrand.size();
    if(IntSize != M_size)
    {
        std::stringstream error;
        error << "Error: dimension of integrand function " << IntSize << " differs from dimension of mesh " << M_size << ".";
        throw std::runtime_error(error.str());
    }
    Real integral = 0;

    for(auto& facet_it : M_mesh.getFractureFacetsIdsVector())
    {
//       Real _volume = M_properties.getProperties(facet_it.getZoneCode()).M_aperture * facet_it.getFacet().area();
        integral += facet_it.getFacet().area()*integrand(facet_it.getIdAsCell());
    }

    return integral;
}

Real Quadrature::integrate(const std::function<Real(Point3D)>& integrand) const
{
    UInt N = M_mesh.getCellsVector().size() + M_mesh.getFractureFacetsIdsVector().size();
    Vector result(N);
    Real integral = 0;
    result = cellIntegrate(integrand);

    for(UInt it = 0; it < N; ++it)
    {
        integral += result (it);
    }

    return integral;
}

Real Quadrature::integrateMatrix(const std::function<Real(Point3D)>& integrand) const
{
    UInt N = M_mesh.getCellsVector().size();
    Vector result(N);
    Real integral = 0;
    result = cellIntegrateMatrix(integrand);

    for(UInt it = 0; it < N; ++it)
    {
        integral += result (it);
    }

    return integral;
}

Real Quadrature::integrateFractures(const std::function<Real(Point3D)>& integrand) const
{
    UInt N = M_mesh.getFractureFacetsIdsVector().size();
    Vector result(N);
    Real integral = 0;
    result = cellIntegrateFractures(integrand);

    for(UInt it = 0; it < N; ++it)
    {
        integral += result (it);
    }

    return integral;
}

Real Quadrature::L2Norm(const Vector& integrand) const
{
    Real integral = 0;
    for(auto& cell_it : M_mesh.getCellsVector())
    {
        integral += cell_it.getVolume()*integrand(cell_it.getId())*integrand(cell_it.getId());
    }

    for(auto& facet_it : M_mesh.getFractureFacetsIdsVector())
    {
//        Real _volume = M_properties.getProperties(facet_it.getZoneCode()).M_aperture * facet_it.getFacet().area();
        integral += facet_it.getFacet().area()*integrand(facet_it.getIdAsCell())*integrand(facet_it.getIdAsCell());
    }

    return sqrt(integral);
}

Real Quadrature::L2NormMatrix(const Vector& integrand) const
{
    Real integral = 0;
    for(auto& cell_it : M_mesh.getCellsVector())
    {
        integral += cell_it.getVolume()*integrand(cell_it.getId())*integrand(cell_it.getId());
    }
    return sqrt(integral);
}

Real Quadrature::L2NormFractures(const Vector& integrand) const
{
    Real integral = 0;

    for(auto& facet_it : M_mesh.getFractureFacetsIdsVector())
    {
//        Real _volume = M_properties.getProperties(facet_it.getZoneCode()).M_aperture * facet_it.getFacet().area();
        integral += facet_it.getFacet().area()*integrand(facet_it.getIdAsCell())*integrand(facet_it.getIdAsCell());
    }

    return sqrt(integral);
}

Vector Quadrature::cellIntegrate (const std::function<Real(Point3D)> & func) const
{
    UInt N = M_mesh.getCellsVector().size() + M_mesh.getFractureFacetsIdsVector().size();
    Vector result(N);

    for(auto& cell_it : M_mesh.getCellsVector())
    {
        result(cell_it.getId()) = M_quadrature->apply(cell_it, func);
    }

    for(auto& facet_it : M_mesh.getFractureFacetsIdsVector())
    {
        result (facet_it.getIdAsCell()) = M_fractureQuadrature->apply(M_mesh.getFacetsVector()[facet_it.getId()], func);
    }

    return result;
}

Vector Quadrature::cellIntegrateMatrix (const std::function<Real(Point3D)> & func) const
{
    UInt N = M_mesh.getCellsVector().size() + M_mesh.getFractureFacetsIdsVector().size();
    Vector result( Vector::Zero(N) );

    for(auto& cell_it : M_mesh.getCellsVector())
    {
        result (cell_it.getId()) = M_quadrature->apply(cell_it, func);
    } // for

    return result;
} // CellIntegrateMatrix

Vector Quadrature::cellIntegrateFractures (const std::function<Real(Point3D)> & func) const
{
    UInt N = M_mesh.getCellsVector().size() + M_mesh.getFractureFacetsIdsVector().size();
    Vector result( Vector::Zero(N) );

    for(auto& facet_it : M_mesh.getFractureFacetsIdsVector())
    {
        result (facet_it.getIdAsCell()) = M_fractureQuadrature->apply(M_mesh.getFacetsVector()[facet_it.getId()], func);
    } // for

    return result;
} // CellIntegrateFractures

} // namespace FVCode3D
