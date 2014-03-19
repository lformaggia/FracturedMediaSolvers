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

template<class Solver, class QRMatrix, class QRFracture>
class DarcyPseudoSteady : public Problem<Solver, QRMatrix, QRFracture>
{
public:

	typedef Geometry::Rigid_Mesh Rigid_Mesh;

	DarcyPseudoSteady(const Rigid_Mesh & mesh, const BoundaryConditions & bc, const Func & f):
		Problem<Solver, QRMatrix, QRFracture>(mesh, bc, f) {};

	virtual void solve();

	virtual ~DarcyPseudoSteady() {};

private:

	DarcyPseudoSteady();

	DarcyPseudoSteady(const DarcyPseudoSteady &);

};

template <class Solver, class QRMatrix, class QRFracture>
void DarcyPseudoSteady< Solver, QRMatrix, QRFracture >::solve()
{

}

#endif /* DARCYPSEUDOSTEADY_HPP_ */
