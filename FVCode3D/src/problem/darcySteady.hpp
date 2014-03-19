/*!
 *	@file darcySteady.hpp
 *	@brief Class that defines and solves the Darcy's problem at steady state.
 */

#ifndef DARCYSTEADY_HPP_
#define DARCYSTEADY_HPP_

#include "problem/problem.hpp"
#include "quadrature/Quadrature.hpp"
#include "solver/solver.hpp"
#include "assembler/stiffness.hpp"

template<class Solver, class QRMatrix, class QRFracture>
class DarcySteady : public Problem<Solver, QRMatrix, QRFracture>
{
public:

	typedef Geometry::Rigid_Mesh Rigid_Mesh;

	DarcySteady(const Rigid_Mesh & mesh, const BoundaryConditions & bc, const Func & f):
		Problem<Solver, QRMatrix, QRFracture>(mesh, bc, f) {};

	virtual void solve();

	virtual ~DarcySteady() {};

private:

	DarcySteady();

	DarcySteady(const DarcySteady &);

};

template <class Solver, class QRMatrix, class QRFracture>
void DarcySteady< Solver, QRMatrix, QRFracture >::solve()
{
	this->M_quadrature.reset( new Quadrature(this->M_mesh, QRMatrix(), QRFracture()) );

	Darcy::StiffMatrix S(this->M_mesh, this->M_bc);
	S.assemble();

	Eigen::VectorXd f(S.getSize());
	f = this->M_quadrature->CellIntegrate(this->M_f);

	SpMat A = S.getMatrix();
	Eigen::VectorXd b;
	b = S.getBCVector() + f;

	this->M_solver.reset( new Solver(A,b) );
	this->M_solver->solve();
}

#endif /* DARCYSTEADY_HPP_ */
