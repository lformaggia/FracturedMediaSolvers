/*!
 *	@file darcyPseudoSteady.hpp
 *	@brief Class that defines and solves the Darcy's problem at pseudo-steady state.
 */

#ifndef DARCYPSEUDOSTEADY_HPP_
#define DARCYPSEUDOSTEADY_HPP_

#include "problem/problem.hpp"
#include "quadrature/Quadrature.hpp"
#include "solver/solver.hpp"
#include "assembler/stiffness.hpp"
#include "assembler/mass.hpp"

class Data;

//! Select the time scheme
/*!
 * @enum TimeScheme
 * It is possible to choose the time scheme: "Implicit" type or "BDF2" type.
 */
enum TimeScheme
{
	Implicit = 0,
	BDF2 = 1
};

//! Class that defines the pseudo-steady-state Darcy problem
/*!
 * @class DarcyPseudoSteady
 * This class defines the pseudo-steady-state Darcy problem given a mesh, boundary conditions, source term and the time parameters.
 * It defines the assemble and the solve method.
 * The first template indicates the Solver used to solve the linear system,
 * the second and third template parameter indicate the quadrature rule for the matrix and fracture respectively,
 * the last template parameter indicates the time scheme to employ.
 *
 * The generic template class is only declared, but not defined.
 * Each time scheme requires a specialization.
 */
template <class Solver, class QRMatrix, class QRFracture, Int TimeScheme>
class DarcyPseudoSteady;

//! Specialization class that defines the pseudo-steady-state Darcy problem for the implicit time scheme
template <class Solver, class QRMatrix, class QRFracture>
class DarcyPseudoSteady< Solver, QRMatrix, QRFracture, Implicit > : public Problem<Solver, QRMatrix, QRFracture>
{
public:

	//! Typedef for Rigid_Mesh
	/*!
	 * @typedef Rigid_Mesh
	 * This type definition permits to treat Geometry::Rigid_Mesh as a Rigid_Mesh
	 */
	typedef Geometry::Rigid_Mesh Rigid_Mesh;

    //! No default constructor
	DarcyPseudoSteady() = delete;

	//! No copy constructor
	DarcyPseudoSteady(const DarcyPseudoSteady &) = delete;

	//! Constructor
	/*!
	 * @param mesh reference to a Geometry::Rigid_mesh
	 * @param bc reference to a BoundaryConditions class
	 * @param func reference to a Func
	 * @param data reference to a Data class
	 */
	DarcyPseudoSteady(const Rigid_Mesh & mesh, const BoundaryConditions & bc, const Func & func, const Data & data):
		Problem<Solver, QRMatrix, QRFracture>(mesh, bc, func, data),
		M_tInit(data.getInitialTime()), M_tEnd(data.getEndTime()), M_tStep(data.getTimeStep()),
		M_x(nullptr), M_M(nullptr), M_nStep(0) {};

	//! Get the previous solution
	/*!
	 * @return constant reference vector to the previous solution
	 */
	const Vector & getOldSolution() const { return M_xOld; }

	//! Initialize the matrices
	/*!
	 * Build the matrices and vectors that are not time dependent, i.e.,
	 * the mass matrix, the stiffness matrix (and BCs vector) and the source vector.
	 * To be called only once.
	 */
	virtual void initialize();

	//! Assemble method
	/*!
	 * Update the previous solution and build the right hand side.
	 * To be called at each time step.
	 */
	virtual void assemble();

	//! Solve method
	/*!
	 * Solve the system at the current time step.
	 * To be called at each time step.
	 * @pre call assemble().
	 */
	virtual void solve();

	//! Destructor
	virtual ~DarcyPseudoSteady() {};

protected:

	//! Initial time
	Real M_tInit;
	//! Final time
	Real M_tEnd;
	//! Time step
	Real M_tStep;

	//! Pointer to a constant vector that contains the current solution
	const Vector * M_x;
	//! Vector that contains the previous solution
	Vector M_xOld;
	//! Vector that contains the source/sink term
	Vector M_f;

	//! Pointer to the mass matrix
	std::unique_ptr<Darcy::MassMatrix> M_M;
	//! Pointer to the stiffness matrix
	std::unique_ptr<Darcy::StiffMatrix> M_S;

	//! Number of time step
	UInt M_nStep;
};

//! Specialization class that defines the pseudo-steady-state Darcy problem for the BDF2 time scheme
template <class Solver, class QRMatrix, class QRFracture>
class DarcyPseudoSteady< Solver, QRMatrix, QRFracture, BDF2 > : public Problem<Solver, QRMatrix, QRFracture>
{
public:

	//! Typedef for Rigid_Mesh
	/*!
	 * @typedef Rigid_Mesh
	 * This type definition permits to treat Geometry::Rigid_Mesh as a Rigid_Mesh
	 */
	typedef Geometry::Rigid_Mesh Rigid_Mesh;

	//! No default constructor
	DarcyPseudoSteady() = delete;

	//! No copy constructor
	DarcyPseudoSteady(const DarcyPseudoSteady &) = delete;

	//! Constructor
	/*!
	 * @param mesh reference to a Geometry::Rigid_mesh
	 * @param bc reference to a BoundaryConditions class
	 * @param func reference to a Func
	 * @param data reference to a Data class
	 */
	DarcyPseudoSteady(const Rigid_Mesh & mesh, const BoundaryConditions & bc, const Func & f, const Data & data):
		Problem<Solver, QRMatrix, QRFracture>(mesh, bc, f, data),
		M_tInit(data.getInitialTime()), M_tEnd(data.getEndTime()), M_tStep(data.getTimeStep()),
		M_x(nullptr), M_M(nullptr), M_nStep(0) {};

	//! Get the previous solution
	/*!
	 * @return constant reference vector to the previous solution
	 */
	const Vector & getOldSolution() const { return M_xOld; }

	//! Get the solution two steps before
	/*!
	 * @return constant reference vector to the solution two steps before
	 */
	const Vector & getOldOldSolution() const { return M_xOldOld; }

	//! Initialize the matrices
	/*!
	 * Build the matrices and vectors that are not time dependent, i.e.,
	 * the mass matrix, the stiffness matrix (and BCs vector) and the source vector.
	 * To be called only once.
	 */
	virtual void initialize();

	//! Assemble method
	/*!
	 * Update the previous solutions and build the right hand side.
	 * To be called at each time step.
	 */
	virtual void assemble();

	//! Solve method
	/*!
	 * Solve the system at the current time step.
	 * To be called at each time step.
	 * @pre call assemble().
	 */
	virtual void solve();

	//! Destructor
	virtual ~DarcyPseudoSteady() {};

protected:

	//! Initial time
	Real M_tInit;
	//! Final time
	Real M_tEnd;
	//! Time step
	Real M_tStep;

	//! Pointer to a constant vector that contains the current solution
	const Vector * M_x;
	//! Vector that contains the previous solution
	Vector M_xOld;
	//! Vector that contains the solution two steps before
	Vector M_xOldOld;
	//! Vector that contains the source/sink term
	Vector M_f;

	//! Pointer to the mass matrix
	std::unique_ptr<Darcy::MassMatrix> M_M;
	//! Pointer to the stiffness matrix
	std::unique_ptr<Darcy::StiffMatrix> M_S;

	//! Number of time step
	UInt M_nStep;
};

template <class Solver, class QRMatrix, class QRFracture>
void DarcyPseudoSteady< Solver, QRMatrix, QRFracture, Implicit >::initialize()
{
	this->M_quadrature.reset( new Quadrature(this->M_mesh, QRMatrix(), QRFracture()) );

	M_M.reset( new Darcy::MassMatrix(this->M_mesh) );
	M_M->assemble();
	M_M->getMatrix() = M_M->getMatrix() / M_tStep;

	M_S.reset( new Darcy::StiffMatrix(this->M_mesh, this->M_bc) );
	M_S->assemble();

	this->M_A = M_S->getMatrix() + M_M->getMatrix();

	M_f.resize(M_S->getSize());
	if (this->M_ssOn == Data::SourceSinkOn::Both || this->M_ssOn == Data::SourceSinkOn::Matrix)
		M_f = this->M_quadrature->CellIntegrateMatrix(this->M_func);
	if (this->M_ssOn == Data::SourceSinkOn::Both || this->M_ssOn == Data::SourceSinkOn::Fractures)
		M_f += this->M_quadrature->CellIntegrateFractures(this->M_func);

	M_xOld = Vector::Constant(M_M->getSize(), 0.);
	M_x = &M_xOld;

	M_nStep = (M_tEnd - M_tInit) / M_tStep;
}

template <class Solver, class QRMatrix, class QRFracture>
void DarcyPseudoSteady< Solver, QRMatrix, QRFracture, Implicit >::assemble()
{
	M_xOld = *M_x;

	this->M_b = M_S->getBCVector() + M_f + M_M->getMatrix() * M_xOld;
}

template <class Solver, class QRMatrix, class QRFracture>
void DarcyPseudoSteady< Solver, QRMatrix, QRFracture, Implicit >::solve()
{
	this->M_solver.reset( new Solver(this->M_A,this->M_b) );
	this->M_solver->solve();

	M_x = &(this->M_solver->getSolution());
}


template <class Solver, class QRMatrix, class QRFracture>
void DarcyPseudoSteady< Solver, QRMatrix, QRFracture, BDF2 >::initialize()
{
	this->M_quadrature.reset( new Quadrature(this->M_mesh, QRMatrix(), QRFracture()) );

	M_M.reset( new Darcy::MassMatrix(this->M_mesh) );
	M_M->assemble();
	M_M->getMatrix() = M_M->getMatrix() / M_tStep;

	M_S.reset( new Darcy::StiffMatrix(this->M_mesh, this->M_bc) );
	M_S->assemble();

	this->M_A = M_S->getMatrix() + (3./2.) * M_M->getMatrix();

	M_f.resize(M_S->getSize());
	if (this->M_ssOn == Data::SourceSinkOn::Both || this->M_ssOn == Data::SourceSinkOn::Matrix)
		M_f = this->M_quadrature->CellIntegrateMatrix(this->M_func);
	if (this->M_ssOn == Data::SourceSinkOn::Both || this->M_ssOn == Data::SourceSinkOn::Fractures)
		M_f += this->M_quadrature->CellIntegrateFractures(this->M_func);

	M_xOldOld = Vector::Constant(M_M->getSize(), 0.);
	M_xOld = Vector::Constant(M_M->getSize(), 0.);
	M_x = &M_xOld;

	M_nStep = (M_tEnd - M_tInit) / M_tStep;
}

template <class Solver, class QRMatrix, class QRFracture>
void DarcyPseudoSteady< Solver, QRMatrix, QRFracture, BDF2 >::assemble()
{
	M_xOldOld = M_xOld;
	M_xOld = *M_x;

	this->M_b = M_S->getBCVector() + M_f + 2 * M_M->getMatrix() * M_xOld - (1./2.) * M_M->getMatrix() * M_xOldOld;
}

template <class Solver, class QRMatrix, class QRFracture>
void DarcyPseudoSteady< Solver, QRMatrix, QRFracture, BDF2 >::solve()
{
	this->M_solver.reset( new Solver(this->M_A,this->M_b) );
	this->M_solver->solve();

	M_x = &(this->M_solver->getSolution());
}

#endif /* DARCYPSEUDOSTEADY_HPP_ */
