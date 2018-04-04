#ifndef FVCODE3D_EVALUATE_HPP
#define FVCODE3D_EVALUATE_HPP

#include <FVCode3D/mesh/RigidMesh.hpp>

namespace FVCode3D
{

Vector evaluateMatrix(const Rigid_Mesh & mesh, const std::function<Real(Point3D)> & func);

Vector evaluateFracture(const Rigid_Mesh & mesh, const std::function<Real(Point3D)> & func);

Vector evaluate(const Rigid_Mesh & mesh, const std::function<Real(Point3D)> & func);

Vector evaluateVelocity(const Rigid_Mesh & mesh, const std::function<Real(Point3D)> & funcx, const std::function<Real(Point3D)> & funcy,
	const std::function<Real(Point3D)> & funcz);

} // namespace FVCode3D

#endif // FVCODE3D_EVALUATE_HPP
