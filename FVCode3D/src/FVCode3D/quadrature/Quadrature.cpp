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

Real Quadrature::integrate(const Vector & integrand) throw()
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
        Real _volume = M_properties.getProperties(facet_it.getZoneCode()).M_aperture * facet_it.getFacet().area();
        integral += _volume*integrand(facet_it.getIdAsCell());
        integral -= _volume/2.*integrand(facet_it.getSeparatedCellsIds()[0]);
        integral -= _volume/2.*integrand(facet_it.getSeparatedCellsIds()[1]);
    }

    return integral;
}

Real Quadrature::integrateMatrix(const Vector & integrand) throw()
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
        Real _volume = M_properties.getProperties(facet_it.getZoneCode()).M_aperture * facet_it.getFacet().area();
        integral -= _volume/2.*integrand(facet_it.getSeparatedCellsIds()[0]);
        integral -= _volume/2.*integrand(facet_it.getSeparatedCellsIds()[1]);
    }

    return integral;
}

Real Quadrature::integrateFractures(const Vector & integrand) throw()
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
        Real _volume = M_properties.getProperties(facet_it.getZoneCode()).M_aperture * facet_it.getFacet().area();
        integral += _volume*integrand(facet_it.getIdAsCell());
    }

    return integral;
}

Real Quadrature::integrate(const std::function<Real(Point3D)>& integrand)
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

Real Quadrature::integrateMatrix(const std::function<Real(Point3D)>& integrand)
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

Real Quadrature::integrateFractures(const std::function<Real(Point3D)>& integrand)
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

Real Quadrature::L2Norm(const Vector& integrand)
{
    Real integral = 0;
    for(auto& cell_it : M_mesh.getCellsVector())
    {
        integral += cell_it.getVolume()*integrand(cell_it.getId())*integrand(cell_it.getId());
    }

    for(auto& facet_it : M_mesh.getFractureFacetsIdsVector())
    {
        Real _volume = M_properties.getProperties(facet_it.getZoneCode()).M_aperture * facet_it.getFacet().area();
        integral += _volume*integrand(facet_it.getIdAsCell())*integrand(facet_it.getIdAsCell());
        integral -= _volume/2.*integrand(facet_it.getSeparatedCellsIds()[0])*integrand(facet_it.getSeparatedCellsIds()[0]);
        integral -= _volume/2.*integrand(facet_it.getSeparatedCellsIds()[1])*integrand(facet_it.getSeparatedCellsIds()[1]);
    }

    return sqrt(integral);
}

Real Quadrature::L2NormMatrix(const Vector& integrand)
{
    Real integral = 0;
    for(auto& cell_it : M_mesh.getCellsVector())
    {
        integral += cell_it.getVolume()*integrand(cell_it.getId())*integrand(cell_it.getId());
    }

    for(auto& facet_it : M_mesh.getFractureFacetsIdsVector())
    {
        Real _volume = M_properties.getProperties(facet_it.getZoneCode()).M_aperture * facet_it.getFacet().area();
        integral -= _volume/2.*integrand(facet_it.getSeparatedCellsIds()[0])*integrand(facet_it.getSeparatedCellsIds()[0]);
        integral -= _volume/2.*integrand(facet_it.getSeparatedCellsIds()[1])*integrand(facet_it.getSeparatedCellsIds()[1]);
    }

    return sqrt(integral);
}

Real Quadrature::L2NormFractures(const Vector& integrand)
{
    Real integral = 0;

    for(auto& facet_it : M_mesh.getFractureFacetsIdsVector())
    {
        Real _volume = M_properties.getProperties(facet_it.getZoneCode()).M_aperture * facet_it.getFacet().area();
        integral += _volume*integrand(facet_it.getIdAsCell())*integrand(facet_it.getIdAsCell());
    }

    return sqrt(integral);
}

Vector Quadrature::cellIntegrate (const std::function<Real(Point3D)> & func)
{
    UInt N = M_mesh.getCellsVector().size() + M_mesh.getFractureFacetsIdsVector().size();
    Vector result(N);
    Real partRes;
    UInt NeighboursId;

    for(auto& cell_it : M_mesh.getCellsVector())
    {
        result(cell_it.getId()) = M_quadrature->apply(cell_it, func);
    }

    for(auto& facet_it : M_mesh.getFractureFacetsIdsVector())
    {
        // contribution from the fracture
        Real _volume = M_properties.getProperties(facet_it.getZoneCode()).M_aperture * facet_it.getFacet().area();

        result (facet_it.getIdAsCell()) = M_fractureQuadrature->apply(M_mesh.getFacetsVector()[facet_it.getId()], _volume, func);

        // remove the contribution of the fracture from the first separated cell
        NeighboursId = M_mesh.getFacetsVector()[facet_it.getId()].getSeparatedCellsIds()[0];

        partRes = M_quadrature->apply(M_mesh.getCellsVector()[NeighboursId], func);

        result (NeighboursId) -= _volume/2. * partRes / M_mesh.getCellsVector()[NeighboursId].getVolume();

        // remove the contribution of the fracture from the second separated cell
        NeighboursId = M_mesh.getFacetsVector()[facet_it.getId()].getSeparatedCellsIds()[1];

        partRes = M_quadrature->apply(M_mesh.getCellsVector()[NeighboursId], func);

        result (NeighboursId) -= _volume/2. * partRes / M_mesh.getCellsVector()[NeighboursId].getVolume();
    }

    return result;
}

Vector Quadrature::cellIntegrateMatrix (const std::function<Real(Point3D)> & func)
{
    UInt N = M_mesh.getCellsVector().size() + M_mesh.getFractureFacetsIdsVector().size();
    Vector result( Vector::Zero(N) );
    Real partRes;
    UInt NeighboursId;

    for(auto& cell_it : M_mesh.getCellsVector())
    {
        result (cell_it.getId()) = M_quadrature->apply(cell_it, func);
    } // for

    for(auto& facet_it : M_mesh.getFractureFacetsIdsVector())
    {
        // contribution from the fracture
        const Real _volume = M_properties.getProperties(facet_it.getZoneCode()).M_aperture * facet_it.getFacet().area();

        // remove the contribution of the fracture from the first separated cell
        NeighboursId = M_mesh.getFacetsVector()[facet_it.getId()].getSeparatedCellsIds()[0];

        partRes = M_quadrature->apply(M_mesh.getCellsVector()[NeighboursId], func);

        result (NeighboursId) -= _volume/2. * partRes / M_mesh.getCellsVector()[NeighboursId].getVolume();

        // remove the contribution of the fracture from the second separated cell
        NeighboursId = M_mesh.getFacetsVector()[facet_it.getId()].getSeparatedCellsIds()[1];

        partRes = M_quadrature->apply(M_mesh.getCellsVector()[NeighboursId], func);

        result (NeighboursId) -= _volume/2. * partRes / M_mesh.getCellsVector()[NeighboursId].getVolume();
    } // for

    return result;

} // CellIntegrateMatrix

Vector Quadrature::cellIntegrateFractures (const std::function<Real(Point3D)> & func)
{
    UInt N = M_mesh.getCellsVector().size() + M_mesh.getFractureFacetsIdsVector().size();
    Vector result( Vector::Zero(N) );

    for(auto& facet_it : M_mesh.getFractureFacetsIdsVector())
    {
        // contribution from the fracture
        const Real _volume = M_properties.getProperties(facet_it.getZoneCode()).M_aperture * facet_it.getFacet().area();

        result (facet_it.getIdAsCell()) = M_fractureQuadrature->apply(M_mesh.getFacetsVector()[facet_it.getId()], _volume, func);
    } // for

    return result;
} // CellIntegrateFractures


} // namespace FVCode3D
