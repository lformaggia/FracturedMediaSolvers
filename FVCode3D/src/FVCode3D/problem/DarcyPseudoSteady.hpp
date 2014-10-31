/*!
 *  @file darcyPseudoSteady.hpp
 *  @brief Class that defines and solves the Darcy's problem at pseudo-steady state.
 */

#ifndef DARCYPSEUDOSTEADY_HPP_
#define DARCYPSEUDOSTEADY_HPP_

#include <FVCode3D/core/Data.hpp>
#include <FVCode3D/problem/Problem.hpp>
#include <FVCode3D/quadrature/Quadrature.hpp>
#include <FVCode3D/assembler/Stiffness.hpp>
#include <FVCode3D/assembler/Mass.hpp>

namespace FVCode3D
{

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
 * The first and the second template parameter indicate the quadrature rule for the matrix and fracture respectively,
 * the last template parameter indicates the time scheme to employ.
 *
 * The generic template class is only declared, but not defined.
 * Each time scheme requires a specialization.
 */
template <class QRMatrix, class QRFracture, typename MatrixType, Int TimeScheme>
class DarcyPseudoSteady;

//! Specialization class that defines the pseudo-steady-state Darcy problem for the implicit time scheme
template <class QRMatrix, class QRFracture, typename MatrixType>
class DarcyPseudoSteady< QRMatrix, QRFracture, MatrixType, Implicit > :
    public Problem< QRMatrix, QRFracture, MatrixType>
{
public:

    //! Typedef for the matrix type
    /*!
     * @typedef Matrix_Type
     */
    typedef MatrixType Matrix_Type;

    //! No default constructor
    DarcyPseudoSteady() = delete;

    //! No copy constructor
    DarcyPseudoSteady(const DarcyPseudoSteady &) = delete;

    //! Constructor
    /*!
     * @param solver string the identify the solver
     * @param mesh reference to a Rigid_mesh
     * @param bc reference to a BoundaryConditions class
     * @param func reference to a Func
     * @param data reference to a Data class
     */
    DarcyPseudoSteady(const std::string solver, const Rigid_Mesh & mesh, const BoundaryConditions & bc,
                      const Func & f, const DataPtr_Type & data):
        Problem<QRMatrix, QRFracture, MatrixType>(solver, mesh, bc, f, data),
        M_tStep(data->getTimeStep()),
        M_x(nullptr), M_isInitialized(false),
        M_M(nullptr), M_S(nullptr) {}

    //! Get the previous solution
    /*!
     * @return constant reference vector to the previous solution
     */
    const Vector & getOldSolution() const { return M_xOld; }

    //! Get the stiffness matrix
    /*!
     * @return the stiffness matrix
     */
    virtual Matrix_Type & getStiffnessMatrix() { return this->M_S->getMatrix(); }

    //! Assemble the linear system
    /*!
     * It calls initialize(), assembleMatrix() and assembleVector()
     */
    virtual void assemble();

    //! Initialize the matrices
    /*!
     * Build the matrices and vectors that are not time dependent, i.e.,
     * the mass matrix, the stiffness matrix, the source vector and the BCs.
     * The BCs and the source are not added to the RHS.
     * To be called only once.
     */
    virtual void initialize();

    //! Assemble matrix method
    /*!
     * Nothing to do
     */
    virtual void assembleMatrix() {};

    //! Assemble vector method
    /*!
     * Update the previous solutions and build the right hand side
     * with BCs, source vector.
     * To be called at each time step.
     */
    virtual void assembleVector();

    //! Solve method
    /*!
     * Solve the system at the current time step.
     * To be called at each time step.
     * @pre call assemble().
     */
    virtual void solve();

    //! Destructor
    virtual ~DarcyPseudoSteady() = default;

protected:

    //! Time step
    Real M_tStep;

    //! Pointer to a constant vector that contains the current solution
    const Vector * M_x;
    //! Vector that contains the previous solution
    Vector M_xOld;
    //! Vector that contains the source/sink term
    Vector M_f;

    //! True if the system is initialized
    bool M_isInitialized;

    //! Pointer to the mass matrix
    std::unique_ptr<MassMatrix> M_M;
    //! Pointer to the stiffness matrix
    std::unique_ptr<StiffMatrix> M_S;
};

//! Specialization class that defines the pseudo-steady-state Darcy problem for the BDF2 time scheme
template <class QRMatrix, class QRFracture, typename MatrixType>
class DarcyPseudoSteady< QRMatrix, QRFracture, MatrixType, BDF2 > :
    public Problem<QRMatrix, QRFracture, MatrixType>
{
public:

    //! Typedef for the matrix type
    /*!
     * @typedef Matrix_Type
     */
    typedef MatrixType Matrix_Type;

    //! No default constructor
    DarcyPseudoSteady() = delete;

    //! No copy constructor
    DarcyPseudoSteady(const DarcyPseudoSteady &) = delete;

    //! Constructor
    /*!
     * @param solver string the identify the solver
     * @param mesh reference to a Rigid_mesh
     * @param bc reference to a BoundaryConditions class
     * @param func reference to a Func
     * @param data reference to a Data class
     */
    DarcyPseudoSteady(const std::string solver, const Rigid_Mesh & mesh, const BoundaryConditions & bc,
                      const Func & f, const DataPtr_Type & data):
        Problem<QRMatrix, QRFracture, MatrixType>(solver, mesh, bc, f, data),
        M_tStep(data->getTimeStep()),
        M_x(nullptr), M_isInitialized(false),
        M_M(nullptr), M_S(nullptr) {};

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

    //! Get the stiffness matrix
    /*!
     * @return the stiffness matrix
     */
    virtual Matrix_Type & getStiffnessMatrix() { return this->M_S->getMatrix(); }

    //! Initialize the matrices
    /*!
     * Build the matrices and vectors that are not time dependent, i.e.,
     * the mass matrix, the stiffness matrix, the source vector and the BCs.
     * The BCs and the source are not added to the RHS.
     * To be called only once.
     */
    virtual void initialize();

    //! Assemble the linear system
    /*!
     * It calls initialize(), assembleMatrix() and assembleVector()
     */
    virtual void assemble();

    //! Assemble matrix method
    /*!
     * Nothing to do
     */
    virtual void assembleMatrix() {};

    //! Assemble vector method
    /*
     * Update the previous solutions and build the right hand side
     * with BCs, source vector.
     * To be called at each time step.
     */
    virtual void assembleVector();

    //! Solve method
    /*!
     * Solve the system at the current time step.
     * To be called at each time step.
     * @pre call assemble().
     */
    virtual void solve();

    //! Destructor
    virtual ~DarcyPseudoSteady() = default;

protected:

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

    //! True if the system is initialized
    bool M_isInitialized;

    //! Pointer to the mass matrix
    std::unique_ptr<MassMatrix> M_M;
    //! Pointer to the stiffness matrix
    std::unique_ptr<StiffMatrix> M_S;
};


//
// Implementation of DarcyPseudoSteady< Implicit >
//

template <class QRMatrix, class QRFracture, typename MatrixType>
void DarcyPseudoSteady< QRMatrix, QRFracture, MatrixType, Implicit >::
assemble()
{
    if(!M_isInitialized)
    {
        initialize();
    }
    assembleMatrix();
    assembleVector();
}

template <class QRMatrix, class QRFracture, typename MatrixType>
void DarcyPseudoSteady< QRMatrix, QRFracture, MatrixType, Implicit >::
initialize()
{
    this->M_quadrature.reset( new Quadrature(this->M_mesh, QRMatrix(), QRFracture()) );

    M_M.reset( new MassMatrix(this->M_mesh) );
    M_M->assemble();
    M_M->closeMatrix();
    M_M->getMatrix() = M_M->getMatrix() / M_tStep;

    M_S.reset( new StiffMatrix(this->M_mesh, this->M_bc) );
    /* MFD */
    M_S->assembleMFD();

    /* FV */
//    M_S->assemble();
//    M_S->closeMatrix();

    this->M_A = M_S->getMatrix() + M_M->getMatrix();

    M_f = Vector::Constant(M_S->getSize(), 0.);
    if (this->M_ssOn == Data::SourceSinkOn::Both || this->M_ssOn == Data::SourceSinkOn::Matrix)
        M_f = this->M_quadrature->cellIntegrateMatrix(this->M_func);
    if (this->M_ssOn == Data::SourceSinkOn::Both || this->M_ssOn == Data::SourceSinkOn::Fractures)
        M_f += this->M_quadrature->cellIntegrateFractures(this->M_func);

    M_xOld = Vector::Constant(M_M->getSize(), 0.);
    M_x = &M_xOld;

    M_isInitialized = true;
} // DarcyPseudoSteady::initialize

template <class QRMatrix, class QRFracture, typename MatrixType>
void DarcyPseudoSteady< QRMatrix, QRFracture, MatrixType, Implicit >::
assembleVector()
{
    M_xOld = *M_x;

    this->M_b = M_S->getBCVector() + M_f + M_M->getMatrix() * M_xOld;
} // DarcyPseudoSteady::assemble

template <class QRMatrix, class QRFracture, typename MatrixType>
void DarcyPseudoSteady< QRMatrix, QRFracture, MatrixType, Implicit >::
solve()
{
    this->M_solver->solve();
    M_x = &(this->M_solver->getSolution());
} // DarcyPseudoSteady::solve


//
// Implementation of DarcyPseudoSteady< BDF2 >
//

template <class QRMatrix, class QRFracture, typename MatrixType>
void DarcyPseudoSteady< QRMatrix, QRFracture, MatrixType, BDF2 >::
assemble()
{
    if(!M_isInitialized)
    {
        initialize();
    }
    assembleMatrix();
    assembleVector();
}

template <class QRMatrix, class QRFracture, typename MatrixType>
void DarcyPseudoSteady< QRMatrix, QRFracture, MatrixType, BDF2 >::
initialize()
{
    this->M_quadrature.reset( new Quadrature(this->M_mesh, QRMatrix(), QRFracture()) );

    M_M.reset( new MassMatrix(this->M_mesh) );
    M_M->assemble();
    M_M->closeMatrix();
    M_M->getMatrix() = M_M->getMatrix() / M_tStep;

    M_S.reset( new StiffMatrix(this->M_mesh, this->M_bc) );

    /* MFD */
    M_S->assembleMFD();

    /* FV */
//    M_S->assemble();
//    M_S->closeMatrix();

    this->M_A = M_S->getMatrix() + (3./2.) * M_M->getMatrix();

    M_f = Vector::Constant(M_S->getSize(), 0.);
    if (this->M_ssOn == Data::SourceSinkOn::Both || this->M_ssOn == Data::SourceSinkOn::Matrix)
        M_f = this->M_quadrature->cellIntegrateMatrix(this->M_func);
    if (this->M_ssOn == Data::SourceSinkOn::Both || this->M_ssOn == Data::SourceSinkOn::Fractures)
        M_f += this->M_quadrature->cellIntegrateFractures(this->M_func);

    M_xOldOld = Vector::Constant(M_M->getSize(), 0.);
    M_xOld = Vector::Constant(M_M->getSize(), 0.);
    M_x = &M_xOld;

    M_isInitialized = true;
} // DarcyPseudoSteady::initialize

template <class QRMatrix, class QRFracture, typename MatrixType>
void DarcyPseudoSteady< QRMatrix, QRFracture, MatrixType, BDF2 >::
assembleVector()
{
    M_xOldOld = M_xOld;
    M_xOld = *M_x;

    this->M_b = M_S->getBCVector() + M_f + 2. * M_M->getMatrix() * M_xOld - (1./2.) * M_M->getMatrix() * M_xOldOld;
} // DarcyPseudoSteady::assemble

template <class QRMatrix, class QRFracture, typename MatrixType>
void DarcyPseudoSteady< QRMatrix, QRFracture, MatrixType, BDF2 >::
solve()
{
    this->M_solver->solve();
    M_x = &(this->M_solver->getSolution());
} // DarcyPseudoSteady::solve

} // namespace FVCode3D
#endif /* DARCYPSEUDOSTEADY_HPP_ */
