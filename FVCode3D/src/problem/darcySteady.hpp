/*!
 *  @file darcySteady.hpp
 *  @brief Class that defines and solves the Darcy's problem at steady state.
 */

#ifndef DARCYSTEADY_HPP_
#define DARCYSTEADY_HPP_

#include "problem/problem.hpp"
#include "quadrature/Quadrature.hpp"
#include "solver/solver.hpp"
#include "assembler/stiffness.hpp"

//! Class that defines the steady-state Darcy problem
/*!
 * @class DarcySteady
 * This class defines the steady-state Darcy problem given a mesh, boundary conditions and source term.
 * It defines the assemble and the solve method.
 * The first template indicates the Solver used to solve the linear system,
 * the second and third template parameter indicate the quadrature rule for the matrix and fracture respectively.
 */
template <class Solver, class QRMatrix, class QRFracture>
class DarcySteady : public Problem<Solver, QRMatrix, QRFracture>
{
public:

    //! Typedef for Rigid_Mesh
    /*!
     * @typedef Rigid_Mesh
     * This type definition permits to treat Geometry::Rigid_Mesh as a Rigid_Mesh
     */
    typedef Geometry::Rigid_Mesh Rigid_Mesh;

    //! No default constructor
    DarcySteady() = delete;

    //! No copy constructor
    DarcySteady(const DarcySteady &) = delete;

    //! Constructor
    /*!
     * @param mesh reference to a Geometry::Rigid_mesh
     * @param bc reference to a BoundaryConditions
     * @param func reference to a Func
     */
    DarcySteady(const Rigid_Mesh & mesh, const BoundaryConditions & bc, const Func & func):
        Problem<Solver, QRMatrix, QRFracture>(mesh, bc, func) {};

    //! Assemble method
    /*!
     * Build the stiffness matrix and the right hand side.
     */
    virtual void assemble();

    //! Solve method
    /*!
     * Solve the system Ax=b.
     * @pre call assemble().
     */
    virtual void solve();

    //! Destructor
    virtual ~DarcySteady() {};
};

template <class Solver, class QRMatrix, class QRFracture>
void DarcySteady< Solver, QRMatrix, QRFracture >::assemble()
{
    this->M_quadrature.reset( new Quadrature(this->M_mesh, QRMatrix(), QRFracture()) );

    Darcy::StiffMatrix S(this->M_mesh, this->M_bc);
    S.assemble();

    Vector f(S.getSize());
    f = this->M_quadrature->CellIntegrate(this->M_func);

    this->M_A = S.getMatrix();
    this->M_b = S.getBCVector() + f;
}

template <class Solver, class QRMatrix, class QRFracture>
void DarcySteady< Solver, QRMatrix, QRFracture >::solve()
{
    this->M_solver.reset( new Solver(this->M_A, this->M_b) );
    this->M_solver->solve();
}

#endif /* DARCYSTEADY_HPP_ */
