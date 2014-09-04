/*!
 * @file BoundingBox.hpp
 * @brief Class that implements a Bounding Box and scaling operations (definitions).
 */

#include <FVCode3D/geometry/BoundingBox.hpp>

namespace FVCode3D
{

BoundingBox::BoundingBox():
    M_xMin(0.), M_xMax(0.), M_yMin(0.), M_yMax(0.), M_zMin(0.), M_zMax(0.), M_L(0.),
    M_mx(1.), M_my(1.), M_mz(1.), M_qx(0.), M_qy(0.), M_qz(0.)
{}

BoundingBox::BoundingBox(const Real xMin, const Real xMax, const Real yMin, const Real yMax, const Real zMin, const Real zMax):
    M_xMin(xMin), M_xMax(xMax), M_yMin(yMin), M_yMax(yMax), M_zMin(zMin), M_zMax(zMax)
{
    computeScalingParameters();
    computeDiagonal();
}

BoundingBox::BoundingBox(const Point3D & pMin, const Point3D & pMax):
    M_xMin(pMin.x()), M_xMax(pMax.x()),
    M_yMin(pMin.y()), M_yMax(pMax.y()),
    M_zMin(pMin.z()), M_zMax(pMax.z())
{
    computeScalingParameters();
    computeDiagonal();
}

BoundingBox::BoundingBox(const std::vector<Point3D> & nodes)
{
    M_xMin = M_yMin = M_zMin = std::numeric_limits<Real>::max();
    M_xMax = M_yMax = M_zMax = std::numeric_limits<Real>::lowest();

    for(auto& p : nodes)
    {
        M_xMin = std::min( M_xMin , p.x() );
        M_yMin = std::min( M_yMin , p.y() );
        M_zMin = std::min( M_zMin , p.z() );

        M_xMax = std::max( M_xMax , p.x() );
        M_yMax = std::max( M_yMax , p.y() );
        M_zMax = std::max( M_zMax , p.z() );
    }

    computeScalingParameters();
    computeDiagonal();
}

void BoundingBox::setExtremes( const Real xmin, const Real xmax,
                               const Real ymin, const Real ymax,
                               const Real zmin, const Real zmax )
{
    M_xMin = xmin;
    M_xMax = xmax;
    M_yMin = ymin;
    M_yMax = ymax;
    M_zMin = zmin;
    M_zMax = zmax;
    computeScalingParameters();
    computeDiagonal();
}

void BoundingBox::computeDiagonal()
{
    M_L =   sqrt(
                (M_xMax-M_xMin) * (M_xMax-M_xMin) +
                (M_yMax-M_yMin) * (M_yMax-M_yMin) +
                (M_zMax-M_zMin) * (M_zMax-M_zMin)
                );
}

void BoundingBox::computeScalingParameters()
{
    M_qx = - M_xMin;
    M_qy = - M_yMin;
    M_qz = - M_zMin;

    M_mx = 1. / ( M_xMax - M_xMin );
    M_my = 1. / ( M_yMax - M_yMin );
    M_mz = 1. / ( M_zMax - M_zMin );

//    M_m = mx;
//    if(my<M_m) M_m=my;
//    if(mz<M_m) M_m=mz;
}

void BoundingBox::scaleNodesToUnit(Point3D & node) const
{
    node.linearTransform( M_mx , M_my , M_mz , M_qx * M_mx , M_qy * M_my , M_qz * M_mz );
}

void BoundingBox::scaleNodesToUnit(std::vector<Point3D> & nodes) const
{
    for(auto& p : nodes)
    {
        p.linearTransform( M_mx , M_my , M_mz , M_qx * M_mx , M_qy * M_my , M_qz * M_mz );
    }
}

void BoundingBox::scaleNodesToPhysical(Point3D & node) const
{
    node.linearTransform( 1. / M_mx , 1. / M_my , 1. / M_mz , - M_qx, - M_qy, - M_qz );
}

void BoundingBox::scaleNodesToPhysical(std::vector<Point3D> & nodes) const
{
    for(auto& p : nodes)
    {
        p.linearTransform( 1. / M_mx , 1. / M_my , 1. / M_mz , - M_qx, - M_qy, - M_qz );
    }
}

} // namespace FVCode3D
