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

//! Class that defines a generic problem
/*!
 *	@class Problem
 *	This class defines a generic problem given a mesh, boundary conditions and source term.
 *	It is an abstract base class. It declares the assemble and the solve method.
 *	The first template indicates the Solver used to solve the linear system,
 *	the second and third template parameter indicate the quadrature rule for the matrix and fracture respectively.
 */
template <class Solver, class QRMatrix, class QRFracture>
class Problem
{
public:

	//! Typedef for Rigid_Mesh
	/*!
	 * @typedef Rigid_Mesh
	 * This type definition permits to treat Geometry::Rigid_Mesh as a Rigid_Mesh.
	 */
	typedef Geometry::Rigid_Mesh Rigid_Mesh;

	//! Constructor
	/*!
	 * @param mesh reference to a Geometry::Rigid_mesh
	 * @param bc reference to a BoundaryConditions
	 * @param func reference to a Func
	 */
	Problem(const Rigid_Mesh & mesh, const BoundaryConditions & bc, const Func & func):
		M_mesh(mesh), M_bc(bc), M_func(func), M_quadrature(nullptr), M_solver(nullptr) {};

	//! @name Get Methods
	//@{
	//! Get the mesh
	/*!
	 * @return a constant reference to the Geometry::Rigid_Mesh
	 */
	const Rigid_Mesh & getMesh() const { return M_mesh; }

	//! Get the boundary conditions
	/*!
	 * @return a constant reference to the boundary conditions
	 */
	const BoundaryConditions & getBC() const { return M_bc; }

	//! Get the source/sink function
	/*!
	 * @return a constant reference to the source/sink function
	 */
	const Func & getF() const { return M_func; }

	//! Get the class Quadrature
	/*!
	 * @return a constant reference to the class Quadrature
	 */
	const Quadrature & getQuadrature() const { return *M_quadrature; }

	//! Get the solver
	/*!
	 * @return a constant reference to the Solver
	 */
	const Solver & getSolver() const { return *M_solver; }
	//@}
	
	//! Assemble method
	virtual void assemble() = 0;
	
	//! Solve method
	virtual void solve() = 0;

	virtual void assembleAndSolve() { assemble(); solve(); }
	
	//! Destructor
	virtual ~Problem() {};

protected:

	//! Constant reference to a Geometry::Rigid_Mesh
	const Rigid_Mesh & M_mesh;
	//! Constant reference to the boundary conditions
	const BoundaryConditions & M_bc;
	//! Constant reference to the source/sink function
	const Func & M_func;
	//! Pointer to the quadrature class
	std::unique_ptr<Quadrature> M_quadrature;
	//! Pointer to the solver class
	std::unique_ptr<Solver> M_solver;

private:
  
	//! No default constructor
	Problem();
	
	//! No copy constructor
	Problem(const Problem &);
};

#endif /* PROBLEM_HPP_ */
