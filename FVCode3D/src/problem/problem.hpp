/*!
 *  @file problem.hpp
 *  @brief Base class that defines a generic problem.
 */

#ifndef PROBLEM_HPP_
#define PROBLEM_HPP_

#include "core/TypeDefinition.hpp"
#include "mesh/Rigid_Mesh.hpp"
#include "boundaryCondition/BC.hpp"

class Quadrature;
class Data;

//! Class that defines a generic problem
/*!
 *  @class Problem
 *  This class defines a generic problem given a mesh, boundary conditions and source term.
 *  It is an abstract base class. It declares the assemble and the solve method.
 *  The first template indicates the Solver used to solve the linear system,
 *  the second and third template parameter indicate the quadrature rule for the matrix and fracture respectively.
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

    //! No default constructor
    Problem() = delete;

    //! No copy constructor
    Problem(const Problem &) = delete;

    //! Constructor
    /*!
     * @param mesh reference to a Geometry::Rigid_mesh
     * @param bc reference to a BoundaryConditions
     * @param func reference to a Func
     */
    Problem(const Rigid_Mesh & mesh, const BoundaryConditions & bc, const Func & func, const Data & data);

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

    //! Get the matrix A
    /*!
     * @return the matrix A
     */
    SpMat & getA() { return M_A; }

    //! Get the RHS
    /*!
     * @return the vector b
     */
    Vector & getb() { return M_b; }
    //@}

    //! Assemble method
    virtual void assemble() = 0;

    //! Solve method
    virtual void solve() = 0;

    //! Assemble and solve
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
    //! Indicates where the source/sink term is applied
    const Data::SourceSinkOn M_ssOn;
    //! Pointer to the quadrature class
    std::unique_ptr<Quadrature> M_quadrature;
    //! Pointer to the solver class
    std::unique_ptr<Solver> M_solver;
    //! Sparse matrix A from the linear system Ax=b
    SpMat M_A;
    //! Vector b from the linear system Ax=b
    Vector M_b;
};

template <class Solver, class QRMatrix, class QRFracture>
Problem< Solver, QRMatrix, QRFracture >::Problem(const Rigid_Mesh & mesh, const BoundaryConditions & bc, const Func & func, const Data & data):
    M_mesh(mesh), M_bc(bc), M_func(func), M_ssOn(data.getSourceSinkOn()), M_quadrature(nullptr), M_solver(nullptr)
{
	M_solver.reset( new Solver() );
}

#endif /* PROBLEM_HPP_ */
