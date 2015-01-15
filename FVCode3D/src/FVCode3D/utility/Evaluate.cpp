#include <FVCode3D/utility/Evaluate.hpp>

namespace FVCode3D
{

Vector evaluateMatrix(const Rigid_Mesh & mesh, const std::function<Real(Point3D)> & func)
{
    const UInt N = mesh.getCellsVector().size();
    Vector result( Vector::Zero(N) );

    for(auto& cell_it : mesh.getCellsVector())
    {
        result(cell_it.getId()) = func(cell_it.getCentroid());
    }

    return result;
} // evaluateMatrix

Vector evaluateFracture(const Rigid_Mesh & mesh, const std::function<Real(Point3D)> & func)
{
    const UInt N = mesh.getFractureFacetsIdsVector().size();
    Vector result( Vector::Zero(N) );

    for(auto& facet_it : mesh.getFractureFacetsIdsVector())
    {
        result(facet_it.getIdAsCell()) = func(facet_it.getCentroid());
    }

    return result;
} // evaluateFracture

Vector evaluate(const Rigid_Mesh & mesh, const std::function<Real(Point3D)> & func)
{
    const UInt N = mesh.getCellsVector().size() + mesh.getFractureFacetsIdsVector().size();
    Vector result( Vector::Zero(N) );

    for(auto& cell_it : mesh.getCellsVector())
    {
        result(cell_it.getId()) = func(cell_it.getCentroid());
    }

    for(auto& facet_it : mesh.getFractureFacetsIdsVector())
    {
        result(facet_it.getIdAsCell()) = func(facet_it.getCentroid());
    }

    return result;
} // evaluate

} // namespace FVCode3D
