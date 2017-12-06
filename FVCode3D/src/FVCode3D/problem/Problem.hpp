/*!
 *  @file problem.hpp
 *  @brief Base class that defines a generic problem.
 */

#ifndef PROBLEM_HPP_
#define PROBLEM_HPP_

#include <utility>

#include <FVCode3D/core/Data.hpp>
#include <FVCode3D/core/TypeDefinition.hpp>
#include <FVCode3D/mesh/RigidMesh.hpp>
#include <FVCode3D/boundaryCondition/BC.hpp>
#include <FVCode3D/solver/SolverHandler.hpp>
#include <FVCode3D/preconditioner/preconditioner.hpp>

namespace FVCode3D
{

class Quadrature;

//! Class that defines a generic problem
/*!
 *  @class Problem
 *  This class defines a generic problem given a mesh, boundary conditions and source term.
 *  It is an abstract base class. It declares the assemble and the solve method.
 *  The first and the second template parameter indicate the quadrature rule for the matrix and fracture respectively.
 */
template <class QRMatrix, class QRFracture>
class Problem
{
public:

//! Typedef for SolverPtr_Type
/*!
 * @typedef SolverPtr_Type
 * This type definition permits to handle a std::tuple<SpMat*,SaddlePointMat*,Vector*> as a AlgebricTuple.
 */
typedef std::tuple<SpMat*,SaddlePointMat*,Vector*> AlgebricTuple;

    //! No default constructor
    Problem() = delete;

    //! No copy constructor
    Problem(const Problem &) = delete;

    //! Constructor
    /*!
     * @param solver string the identify the solver
     * @param mesh reference to a Rigid_mesh
     * @param bc reference to a BoundaryConditions
     * @param func reference to a Func
     * @param data reference to a Data class
     */
    Problem(const std::string solver, const Rigid_Mesh & mesh, const BoundaryConditions & bc,
            const Func & func, const DataPtr_Type & data);

    //! @name Get Methods
    //@{
    //! Get the mesh
    /*!
     * @return a constant reference to the Rigid_Mesh
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

    //! Get where the source/sink term is applied
    /*!
     * @return where the source/sink term is applied
     */
    Data::SourceSinkOn getSourceSinkOn() const { return M_ssOn; }

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

    //! Get the solver
    /*!
     * @return a reference to the Solver
     */
    Solver & getSolver() { return *M_solver; }

    //! Get the solver pointer
    /*!
     * @return a pointer to the Solver
     */
    Solver * getSolverPtr() { return M_solver.get(); }
    
    //! Get the algebraic counterpart of the problem
    /*!
     * @return a tuple of three different pointers.
     * The first pointer is to an SpMat that, if the solver is of direct type,
     * is the system matrix, otherwise it is a null ptr.
     * The second pointer is to an SaddlePointMat that, if the solver is of iterative type, 
     * is the system block matrix, otherwise it is a null ptr.
     * The third ptr is to a Vector that is the rhs of the system.
     */
    AlgebricTuple getAlgebry() 
    { 
		SpMat          *   M_A(nullptr);
		SaddlePointMat *   M_SP(nullptr);
		Vector         *   M_b(& this->M_solver->getb());
    
		if(dynamic_cast<DirectSolver*>(this->getSolverPtr()))
			M_A = & dynamic_cast<DirectSolver*>(this->getSolverPtr())->getA();
        
		if(dynamic_cast<IterativeSolver*>(this->getSolverPtr()))
			M_SP = & dynamic_cast<IterativeSolver*>(this->getSolverPtr())->getA();   
			
		return std::make_tuple(M_A,M_SP,M_b);
	}
    //@}

    //! Assemble method
    /*!
     * It calls assembleMatrix() and assembleVector()
     */
    void assemble() { assembleMatrix(); assembleVector(); };

    //! Assemble matrix method
    virtual void assembleMatrix() = 0;

    //! Assemble vector method
    virtual void assembleVector() = 0;

    //! Solve method
    virtual void solve() = 0;

    //! Assemble and solve
    /*!
     * It calls assemble() and solve()
     */
    virtual void assembleAndSolve() { assemble(); solve(); }

    //! Destructor
    virtual ~Problem() = default;

protected:

    //! Constant reference to a Rigid_Mesh
    const Rigid_Mesh & M_mesh;
    //! Constant reference to the boundary conditions
    const BoundaryConditions & M_bc;
    //! Constant reference to the source/sink function
    const Func & M_func;
    //! Indicates where the source/sink term is applied
    const Data::SourceSinkOn M_ssOn;
    //! Indicates the numerical method
    const Data::NumericalMethodType M_numet;
    //! Pointer to the quadrature class
    std::unique_ptr<Quadrature> M_quadrature;
    //! Pointer to the solver class
    SolverPtr_Type M_solver;
};

template <class QRMatrix, class QRFracture>
Problem<QRMatrix, QRFracture>::
Problem(const std::string solver, const Rigid_Mesh & mesh, const BoundaryConditions & bc,
        const Func & func, const DataPtr_Type & data):
    M_mesh(mesh),
    M_bc(bc),
    M_func(func),
    M_ssOn(data->getSourceSinkOn()),
    M_numet(data->getNumericalMethodType()),
    M_quadrature(nullptr),
    M_solver( SolverHandler::Instance().getProduct(solver) )
	{} // Problem

} // namespace FVCode3D

#endif /* PROBLEM_HPP_ */
