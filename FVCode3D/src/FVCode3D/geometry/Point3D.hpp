/*!
 * @file Point3D.hpp
 * @brief Class for a 3D Point.
 */

#ifndef POINT3D_HPP_
#define POINT3D_HPP_

#include <FVCode3D/core/BasicType.hpp>

namespace FVCode3D
{

class CoordinateSystem3D;

//! Class that defines a 3D Point and some operations.
/*!
 * @class Point3D
 * This class defines a point as three Real coordinates.
 * It also implements some operation between two points or between a point and a scalar value.
 */
class Point3D
{
public:

    //! Empty constructor
    /*!
     * Initialize with null coordinates
     */
    Point3D(): M_x(0.) ,M_y(0.) , M_z(0.) {};

    //! Constructor
    /*!
     * Initialize with coordinates
     * @param x x-coordinate
     * @param y y-coordinate
     * @param z z-coordinate
     */
    Point3D(Real x, Real y, Real z): M_x(x) ,M_y(y) , M_z(z) {};

    //! Copy constructor
    /*!
     * @param p point
     */
    Point3D(const Point3D & p): M_x(p.x()), M_y(p.y()), M_z(p.z()) {};

    //! Get x-coordinate (const)
    /*!
     * @return x-coordinate
     */
    Real x() const { return M_x; };

    //! Get y-coordinate (const)
    /*!
     * @return x-coordinate
     */
    Real y() const { return M_y; };

    //! Get x-coordinate (const)
    /*!
     * @return x-coordinate
     */
    Real z() const { return M_z; };

    //! Get reference to x-coordinate
    /*!
     * @return x-coordinate
     */
    Real & x() { return M_x; };

    //! Get reference to y-coordinate
    /*!
     * @return y-coordinate
     */
    Real & y() { return M_y; };

    //! Get reference to z-coordinate
    /*!
     * @return z-coordinate
     */
    Real & z() { return M_z; };

    //! Compute norm of the current point
    /*!
     * @return norm of the point
     */
    Real norm() const { return std::sqrt( M_x*M_x + M_y*M_y + M_z*M_z );}

    //! Normalize the current point
    void normalize();

    //! Return the normalized point
    Point3D normalized() const { return Point3D( *this ) / norm(); }

    //! Set coordinates
    /*!
     * @param x x-coordinate
     * @param y y-coordinate
     * @param z z-coordinate
     */
    void setValues(const Real x, const Real y, const Real z);

    //! Apply a linear transformation to the point
    /*!
     * @param scaling scale factor
     * @param xShift translation along x
     * @param yShift translation along y
     * @param zShift translation along z
     */
    void linearTransform(const Real scaling, const Real xShift, const Real yShift, const Real zShift);

    //! Apply a linear transformation to the point
    /*!
     * @param xScaling scale factor along x
     * @param yScaling scale factor along y
     * @param zScaling scale factor along z
     * @param xShift translation along x
     * @param yShift translation along y
     * @param zShift translation along z
     */
    void linearTransform(const Real xScaling, const Real yScaling, const Real zScaling,
                         const Real xShift, const Real yShift, const Real zShift);

    //! Convert the point to a different coordinate system
    /*!
     * @param coordSys new coordinate system
     * @param origin origin of the coordinate system
     * @return the point in the new coordinate system
     */
    Point3D convertInLocalCoordinates(const CoordinateSystem3D & coordSys, const Point3D & origin) const;

    //! Assignment operator
    Point3D & operator=(const Point3D & p);

    //! += operator
    Point3D & operator+=(const Point3D & p);

    //! -= operator
    Point3D & operator-=(const Point3D & p);

    //! *= operator
    Point3D & operator*=(const Real & r);

    //! /= operator
    Point3D & operator/=(const Real & r);

    //! Access to the coord-th component of the point (const)
    /*!
     * @param coord subscript of the coordinate (0, 1 ,2)
     * @return the coord-th component
     */
    Real operator[](const UInt coord) const throw();

    //! Access to the coord-th component of the point
    /*!
     * @param coord subscript of the coordinate (0, 1 ,2)
     * @return the coord-th component
     */
    Real & operator[](const UInt coord) throw();

    //! Get tolerance on x
    /*!
     * @return tolerance on x
     */
    static const Real & getToleranceX()
        { return Point3D::S_toleranceX; }

    //! Get tolerance on y
    /*!
     * @return tolerance on y
     */
    static const Real & getToleranceY()
        { return Point3D::S_toleranceY; }

    //! Get tolerance on z
    /*!
     * @return tolerance on z
     */
    static const Real & getToleranceZ()
        { return Point3D::S_toleranceZ; }

    //! Set tolerance on x
    /*!
     * @param tolerance tolerance on x
     */
    static void setToleranceX(const Real tolerance)
        { Point3D::S_toleranceX = tolerance; }

     //! Set tolerance on y
     /*!
      * @param tolerance tolerance on y
      */
    static void setToleranceY(const Real tolerance)
        { Point3D::S_toleranceY = tolerance; }

     //! Set tolerance on z
     /*!
      * @param tolerance tolerance on z
      */
    static void setToleranceZ(const Real tolerance)
        { Point3D::S_toleranceZ = tolerance; }

    //! + operator
    friend Point3D operator+(const Point3D & p1, const Point3D & p2);

    //! - operator
    friend Point3D operator-(const Point3D & p1, const Point3D & p2);

    //! Change sign operator
    friend Point3D operator-(const Point3D & p);

    //! * operator between points (dot product)
    friend Real operator*(const Point3D & p1, const Point3D & p2);

    //! * operator between point and scalar
    friend Point3D operator*(const Point3D & p, const Real r);

    //! * operator between scalar and point
    friend Point3D operator*(const Real r, const Point3D & p);

    //! / operator between point and scalar
    friend Point3D operator/(const Point3D & p, const Real r);

    //! Dot product
    friend Real dotProduct(const Point3D & p1, const Point3D & p2);

    //! Compute inner angle between two vector in radians
    friend Real innerAngleRad(const Point3D & p1, const Point3D & p2);

    //! Compute inner angle between two vector in degrees
    friend Real innerAngleDeg(const Point3D & p1, const Point3D & p2);

    //! Compute distance between points
    friend Real distance(const Point3D & p1, const Point3D & p2);

    //! Cross product
    friend Point3D crossProduct(const Point3D & p1, const Point3D & p2);

    //! Print information about the point
    friend std::ostream & operator<<(std::ostream & os, const Point3D & p);

    //! Relative x-tolerance
    static Real S_toleranceX;
    //! Relative y-tolerance
    static Real S_toleranceY;
    //! Relative z-tolerance
    static Real S_toleranceZ;

private:

    //! x-coordinate
    Real M_x;
    //! y-coordinate
    Real M_y;
    //! z-coordinate
    Real M_z;
};

//! Operator less between points
bool operator<(const Point3D & p1, const Point3D & p2);

} // namespace FVCode3D

#endif /* POINT3D_HPP_ */
