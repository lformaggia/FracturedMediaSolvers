/*!
 *  @file darcySteady.hpp
 *  @brief Class that defines and solves the Darcy's problem at steady state.
 */

#ifndef DARCYSTEADY_HPP_
#define DARCYSTEADY_HPP_

#include "core/Data.hpp"
#include "problem/Problem.hpp"
#include "quadrature/Quadrature.hpp"
#include "solver/Solver.hpp"
#include "assembler/Stiffness.hpp"

namespace FVCode3D
{

//! Class that defines the steady-state Darcy problem
/*!
 * @class DarcySteady
 * This class defines the steady-state Darcy problem given a mesh, boundary conditions and source term.
 * It defines the assemble and the solve method.
 * The first template indicates the Solver used to solve the linear system,
 * the second and third template parameter indicate the quadrature rule for the matrix and fracture respectively.
 */
template <class Solver, class QRMatrix, class QRFracture, typename MatrixType = SpMat>
class DarcySteady : public Problem<Solver, QRMatrix, QRFracture, MatrixType>
{
public:

    //! Typedef for the matrix type
    /*!
     * @typedef Matrix_Type
     */
    typedef MatrixType Matrix_Type;

    //! No default constructor
    DarcySteady() = delete;

    //! No copy constructor
    DarcySteady(const DarcySteady &) = delete;

    //! Constructor
    /*!
     * @param mesh reference to a Rigid_mesh
     * @param bc reference to a BoundaryConditions
     * @param func reference to a Func
     */
    DarcySteady(const Rigid_Mesh & mesh, const BoundaryConditions & bc, const Func & func, const DataPtr_Type & data):
        Problem<Solver, QRMatrix, QRFracture, MatrixType>(mesh, bc, func, data) {};

    //! Assemble matrix method
    /*!
     * Build the stiffness matrix.
     */
    virtual void assembleMatrix();

    //! Assemble vector method
    /*!
     * Build the right hand side.
     */
    virtual void assembleVector();

    //! Solve method
    /*!
     * Solve the system Ax=b.
     * @pre call assemble().
     */
    virtual void solve() { this->M_solver->solve(); }

    //! Destructor
    virtual ~DarcySteady() = default;
};

template <class Solver, class QRMatrix, class QRFracture, typename MatrixType>
void DarcySteady< Solver, QRMatrix, QRFracture, MatrixType >::assembleMatrix()
{
    this->M_quadrature.reset( new Quadrature(this->M_mesh, QRMatrix(), QRFracture()) );

    StiffMatrix S(this->M_mesh, this->M_bc);
    S.assemble();

    this->M_A = S.getMatrix();
    this->M_b = S.getBCVector();
} // DarcySteady::assembleMatrix

template <class Solver, class QRMatrix, class QRFracture, typename MatrixType>
void DarcySteady< Solver, QRMatrix, QRFracture, MatrixType >::assembleVector()
{
    this->M_quadrature.reset( new Quadrature(this->M_mesh, QRMatrix(), QRFracture()) );

    Vector f( Vector::Constant( this->M_A.rows(), 0.) );
    if ( this->M_mesh.getCellsVector().size() != 0
            &&
            ( this->M_ssOn == Data::SourceSinkOn::Both
                ||
              this->M_ssOn == Data::SourceSinkOn::Matrix) )
    {
        f = this->M_quadrature->cellIntegrateMatrix(this->M_func);
    } // if

    if ( this->M_mesh.getFractureFacetsIdsVector().size() != 0
            &&
            ( this->M_ssOn == Data::SourceSinkOn::Both
                ||
              this->M_ssOn == Data::SourceSinkOn::Fractures ) )
    {
        f = this->M_quadrature->cellIntegrateFractures(this->M_func);
    } // if

    if ( this->M_b.size() == 0 )
    {
        this->M_b = f;
    } // if
    else
    {
        this->M_b += f;
    } // else
} // DarcySteady::assembleVector

} // namespace FVCode3D
#endif /* DARCYSTEADY_HPP_ */