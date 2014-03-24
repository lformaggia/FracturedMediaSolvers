/*!
 *	@file problem.hpp
 *	@brief Base class that defines a generic problem.
 */

#ifndef PROBLEM_HPP_
#define PROBLEM_HPP_

#include "core/TypeDefinition.hpp"
#include "mesh/Rigid_Mesh.hpp"
#include "boundaryCondition/BC.hpp"

class Quadrature;

template <class Solver, class QRMatrix, class QRFracture>
class Problem
{
public:

	typedef Geometry::Rigid_Mesh Rigid_Mesh;

	Problem(const Rigid_Mesh & mesh, const BoundaryConditions & bc, const Func & func):
		M_mesh(mesh), M_bc(bc), M_func(func), M_quadrature(nullptr), M_solver(nullptr) {};

	const Rigid_Mesh & getMesh() const { return M_mesh; }

	const BoundaryConditions & getBC() const { return M_bc; }

	const Func & getF() const { return M_func; }

	const Quadrature & getQuadrature() const { return *M_quadrature; }

	const Solver & getSolver() const { return *M_solver; }

	virtual void solve() = 0;

	virtual ~Problem() {};

protected:

	const Rigid_Mesh & M_mesh;
	const BoundaryConditions & M_bc;
	const Func & M_func;
	std::unique_ptr<Quadrature> M_quadrature;
	std::unique_ptr<Solver> M_solver;

private:
	Problem();
	Problem(const Problem &);

};

#endif /* PROBLEM_HPP_ */
