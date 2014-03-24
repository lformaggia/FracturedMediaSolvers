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

enum TimeScheme
{
	Implicit,
	BDF2
};

template <class Solver, class QRMatrix, class QRFracture, Int TimeScheme>
class DarcyPseudoSteady;

template <class Solver, class QRMatrix, class QRFracture>
class DarcyPseudoSteady< Solver, QRMatrix, QRFracture, Implicit > : public Problem<Solver, QRMatrix, QRFracture>
{
public:

	typedef Geometry::Rigid_Mesh Rigid_Mesh;

	DarcyPseudoSteady(const Rigid_Mesh & mesh, const BoundaryConditions & bc, const Func & func, const Data & data):
		Problem<Solver, QRMatrix, QRFracture>(mesh, bc, func),
		M_tInit(data.getInitialTime()), M_tEnd(data.getEndTime()), M_tStep(data.getTimeStep()),
		M_x(nullptr), M_M(nullptr), M_nStep(0) {};

	const Vector & getOldSolution() const { return M_xOld; }

	virtual void initialize();

	virtual void solve();

	virtual ~DarcyPseudoSteady() {};

protected:

	Real M_tInit;
	Real M_tEnd;
	Real M_tStep;

	const Vector * M_x;
	Vector M_xOld;
	Vector M_f;
	SpMat M_A;
	std::unique_ptr<Darcy::MassMatrix> M_M;
	std::unique_ptr<Darcy::StiffMatrix> M_S;
	UInt M_nStep;

private:

	DarcyPseudoSteady();

	DarcyPseudoSteady(const DarcyPseudoSteady &);

};

template <class Solver, class QRMatrix, class QRFracture>
class DarcyPseudoSteady< Solver, QRMatrix, QRFracture, BDF2 > : public Problem<Solver, QRMatrix, QRFracture>
{
public:

	typedef Geometry::Rigid_Mesh Rigid_Mesh;

	DarcyPseudoSteady(const Rigid_Mesh & mesh, const BoundaryConditions & bc, const Func & f, const Data & data):
		Problem<Solver, QRMatrix, QRFracture>(mesh, bc, f),
		M_tInit(data.getInitialTime()), M_tEnd(data.getEndTime()), M_tStep(data.getTimeStep()),
		M_x(nullptr), M_M(nullptr), M_nStep(0) {};

	const Vector & getOldSolution() const { return M_xOld; }

	const Vector & getOldOldSolution() const { return M_xOldOld; }

	virtual void initialize();

	virtual void solve();

	virtual ~DarcyPseudoSteady() {};

protected:

	Real M_tInit;
	Real M_tEnd;
	Real M_tStep;

	const Vector * M_x;
	Vector M_xOld;
	Vector M_xOldOld;
	Vector M_f;
	SpMat M_A;
	std::unique_ptr<Darcy::MassMatrix> M_M;
	std::unique_ptr<Darcy::StiffMatrix> M_S;
	UInt M_nStep;

private:

	DarcyPseudoSteady();

	DarcyPseudoSteady(const DarcyPseudoSteady &);

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

	M_A = M_S->getMatrix() + M_M->getMatrix();

	M_f.resize(M_S->getSize());
	M_f = this->M_quadrature->CellIntegrate(this->M_func);

	M_xOld = Vector::Constant(M_M->getSize(), 0.);
	M_x = &M_xOld;

	M_nStep = (M_tEnd - M_tInit) / M_tStep;
}

template <class Solver, class QRMatrix, class QRFracture>
void DarcyPseudoSteady< Solver, QRMatrix, QRFracture, Implicit >::solve()
{
	M_xOld = *M_x;

	Vector b(M_M->getSize());
	b = M_S->getBCVector() + M_f + M_M->getMatrix() * M_xOld;

	this->M_solver.reset( new Solver(M_A,b) );
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

	M_A = M_S->getMatrix() + (3./2.) * M_M->getMatrix();

	std::cout << std::endl;
	std::cout<<"S norm: " << M_S->getMatrix().norm() << std::endl;
	std::cout<<"A norm: " << M_A.norm() << std::endl;

	M_f.resize(M_S->getSize());
	M_f = this->M_quadrature->CellIntegrate(this->M_func);

	std::cout<<"func norm: " << M_f.norm() << std::endl;

	M_xOldOld = Vector::Constant(M_M->getSize(), 0.);
	M_xOld = Vector::Constant(M_M->getSize(), 0.);
	M_x = &M_xOld;

	std::cout<<"Sol OldOld norm: " << M_xOldOld.norm() << std::endl;
	std::cout<<"Sol Old norm: " << M_xOld.norm() << std::endl;
	std::cout<<"Sol norm: " << M_x->norm() << std::endl;

	M_nStep = (M_tEnd - M_tInit) / M_tStep;
}

template <class Solver, class QRMatrix, class QRFracture>
void DarcyPseudoSteady< Solver, QRMatrix, QRFracture, BDF2 >::solve()
{
	M_xOldOld = M_xOld;
	M_xOld = *M_x;

	Vector b(M_M->getSize());
	b = M_S->getBCVector() + M_f + 2 * M_M->getMatrix() * M_xOld - (1./2.) * M_M->getMatrix() * M_xOldOld;

	this->M_solver.reset( new Solver(M_A,b) );
	this->M_solver->solve();

	M_x = &(this->M_solver->getSolution());

	std::cout<<"S norm: " << M_S->getMatrix().norm() << std::endl;
	std::cout<<"A norm: " << M_A.norm() << std::endl;

	std::cout<<"Sol OldOld norm: " << M_xOldOld.norm() << std::endl;
	std::cout<<"Sol Old norm: " << M_xOld.norm() << std::endl;
	std::cout<<"Sol norm: " << M_x->norm() << std::endl;
}

#endif /* DARCYPSEUDOSTEADY_HPP_ */
