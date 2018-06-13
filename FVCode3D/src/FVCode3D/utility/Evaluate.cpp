#include <FVCode3D/utility/Evaluate.hpp>
#include <cmath>
#include <Eigen/LU>
#include <unsupported/Eigen/SparseExtra>

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
    const UInt M = mesh.getCellsVector().size();
    Vector result( Vector::Zero(N) );

    for(auto& facet_it : mesh.getFractureFacetsIdsVector())
    {
        result(facet_it.getIdAsCell() - M) = func(facet_it.getCentroid());
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

Vector evaluateVelocity(const Rigid_Mesh & mesh, const std::function<Real(Point3D)> & funcx, const std::function<Real(Point3D)> & funcy,
	const std::function<Real(Point3D)> & funcz)
{
	
UInt numFacetsTot = mesh.getFacetsVector().size() + mesh.getFractureFacetsIdsVector().size();
Vector result(Vector::Zero(numFacetsTot));

for(auto & facet : mesh.getFacetsVector())
{
	auto normal = facet.getUnsignedNormal();
	auto & M_vertexIds = facet.getVerticesIds();
	assert(M_vertexIds.size() >= 3);
	
    if(M_vertexIds.size() > 3)
    {
        Point3D center(0., 0., 0.);
        for(auto id : M_vertexIds)
        {
            center += facet.getMesh()->getNodesVector()[id];
        }
        center /= M_vertexIds.size();

        auto it1 = M_vertexIds.begin();
        auto it2 = it1;
        it2++;

        for( ; it2 != M_vertexIds.end(); ++it1, ++it2)
        {
            const Point3D B = facet.getMesh()->getNodesVector()[*it1];
            const Point3D C = facet.getMesh()->getNodesVector()[*it2];

//			Point3D normTr = FVCode3D::computeNormal(center, B, C);
            Real area = FVCode3D::triangleArea(center, B, C);
            Point3D TriaCentroid = FVCode3D::triangleCentroid(center,B,C);
            result[facet.getId()] += area * ( funcx(TriaCentroid)*normal.x() + funcy(TriaCentroid)*normal.y() + funcz(TriaCentroid)*normal.z() );
        }

        const Point3D B = facet.getMesh()->getNodesVector()[ M_vertexIds[M_vertexIds.size()-1] ];
        const Point3D C = facet.getMesh()->getNodesVector()[ M_vertexIds[0] ];

//		Point3D normTr = FVCode3D::computeNormal(center, B, C);
        Real area = FVCode3D::triangleArea(center, B, C);
        Point3D TriaCentroid = FVCode3D::triangleCentroid(center,B,C);
        result[facet.getId()] += area * ( funcx(TriaCentroid)*normal.x() + funcy(TriaCentroid)*normal.y() + funcz(TriaCentroid)*normal.z() );

        result[facet.getId()] /= facet.area();
    }
    else
    {
        const Point3D A(facet.getMesh()->getNodesVector()[M_vertexIds[0]]);
        const Point3D B(facet.getMesh()->getNodesVector()[M_vertexIds[1]]);
        const Point3D C(facet.getMesh()->getNodesVector()[M_vertexIds[2]]);
//		Point3D normTr = FVCode3D::computeNormal(A, B, C);
        Real area = FVCode3D::triangleArea(A, B, C);
        Point3D TriaCentroid = FVCode3D::triangleCentroid(A, B, C);
        result[facet.getId()] = funcx(TriaCentroid)*normal.x() + funcy(TriaCentroid)*normal.y() + funcz(TriaCentroid)*normal.z();
        
        result[facet.getId()] /= facet.area();	
    }	
}
return result;	
}

} // namespace FVCode3D
