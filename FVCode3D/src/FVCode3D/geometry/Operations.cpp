/*!
 *  @file operations.cpp
 *  @brief Some useful operations (definitions).
 */

#include <FVCode3D/geometry/Operations.hpp>
#include <FVCode3D/geometry/Point3D.hpp>

namespace FVCode3D
{

Point3D computeNormal( const Point3D & A, const Point3D & B, const Point3D & C)
{
    Point3D n( crossProduct(A-C,B-C) );
    n.normalize();
    return n;
}

Point3D rotateOf(const Point3D & p, const Real angleDeg)
{
    Real angleRad = angleDeg / 180 * _PI_;
    Point3D r1(std::cos(angleRad), std::sin(angleRad), 0.);
    Point3D r2(-std::sin(angleRad), std::cos(angleRad), 0.);

    Real x, y;
    x = dotProduct(r1,p);
    y = dotProduct(r2,p);

    return Point3D(x,y,p.z());
}

Real triangleArea(const Point3D & A, const Point3D & B, const Point3D & C)
{
    return crossProduct(B-A,C-A).norm() / 2.;
}

Point3D triangleCentroid(const Point3D & A, const Point3D & B, const Point3D & C)
{
    return (A + B + C) / 3.;
}

Real tetrahedronVolume(const std::vector<Point3D> & nodes)
{
    return std::fabs( dotProduct( nodes[2]-nodes[3] , crossProduct(nodes[0]-nodes[3],nodes[1]-nodes[3]) ) ) / 6.;
}

void computeBoundingBox(const std::vector<Point3D> & nodes, Point3D & pMin, Point3D & pMax)
{
    pMin.setValues(std::numeric_limits<Real>::max(),std::numeric_limits<Real>::max(),std::numeric_limits<Real>::max());
    pMax.setValues(std::numeric_limits<Real>::lowest(),std::numeric_limits<Real>::lowest(),std::numeric_limits<Real>::lowest());

    Real & xMin = pMin.x();
    Real & yMin = pMin.y();
    Real & zMin = pMin.z();
    Real & xMax = pMax.x();
    Real & yMax = pMax.y();
    Real & zMax = pMax.z();

    for(auto p : nodes)
    {
        xMin = std::min( xMin , p.x() );
        yMin = std::min( yMin , p.y() );
        zMin = std::min( zMin , p.z() );

        xMax = std::max( xMax , p.x() );
        yMax = std::max( yMax , p.y() );
        zMax = std::max( zMax , p.z() );
    }
}

}// namespace FVCode3D
