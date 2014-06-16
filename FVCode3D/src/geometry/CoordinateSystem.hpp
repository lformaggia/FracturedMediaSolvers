/*!
 * @file CoordinateSystem.hpp
 * @brief Class for a 3D Cartesian coordinate system.
 */

#ifndef COORDINATESYSTEM_HPP_
#define COORDINATESYSTEM_HPP_

#include "core/TypeDefinition.hpp"
#include "geometry/Point3D.hpp"

namespace FVCode3D
{

class Point3D;

//! Class that implements a Cartesian Coordinate System
/*!
 * @class CoordinateSystem3D
 * This class implements a Cartesian Coordinate System in a 3D space by fixing three directions.
 */
class CoordinateSystem3D{
public:

	//! Empty constructor
	CoordinateSystem3D(): M_u(), M_v(), M_w() {};

	//! Constructor
	/*!
	 * @param u point that sets the x-direction
	 * @param v point that sets the y-direction
	 * @param w point that sets the z-direction
	 */
	CoordinateSystem3D(const Point3D & u, const Point3D & v, const Point3D & w): M_u(u), M_v(v), M_w(w)
				{ M_u.normalize(); M_v.normalize(); M_w.normalize(); };

	//! Get x-direction
	/*!
	 * @return the point that represent the x-direction
	 */
	Point3D getU() const { return M_u; };

	//! Get y-direction
	/*!
	 * @return the point that represent the y-direction
	 */
	Point3D getV() const { return M_v; };

	//! Get z-direction
	/*!
	 * @return the point that represent the z-direction
	 */
	Point3D getW() const { return M_w; };

	//! Compute the Coordinate System moving from a normal that represents the z-direction
	/*!
	 * @param z normal that represent the z-direction
	 */
	void computeCartesianCoordinateSystem(const Point3D & z);

private:

	//! x-direction
	Point3D M_u;
	//! y-direction
	Point3D M_v;
	//! z-direction
	Point3D M_w;
};

} // namespace FVCode3D

#endif /* COORDINATESYSTEM_HPP_ */
