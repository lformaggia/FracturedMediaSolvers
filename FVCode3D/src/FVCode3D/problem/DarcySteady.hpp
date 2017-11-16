/*!
 *  @file DarcySteady.hpp
 *  @brief Class that defines and solves the Darcy's problem at steady state.
 */

#ifndef DARCYSTEADY_HPP_
#define DARCYSTEADY_HPP_

#include <FVCode3D/core/Data.hpp>
#include <FVCode3D/problem/Problem.hpp>
#include <FVCode3D/quadrature/Quadrature.hpp>
#include <FVCode3D/assembler/Stiffness.hpp>
#include <unsupported/Eigen/SparseExtra>

namespace FVCode3D
{

//! Class that defines the steady-state Darcy problem
/*!
 * @class DarcySteady
 * This class defines the steady-state Darcy problem given a mesh, boundary conditions and source term.
 * It defines the assemble and the solve method.
 * The first and the second template parameter indicate the quadrature rule for the matrix and fracture respectively.
 */
template <class QRMatrix, class QRFracture, typename MatrixType = SpMat>
class DarcySteady : public Problem<QRMatrix, QRFracture, MatrixType>
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
     * @param solver string the identify the solver
     * @param mesh reference to a Rigid_mesh
     * @param bc reference to a BoundaryConditions
     * @param func reference to a Func
     * @param data reference to a Data class
     */
    DarcySteady(const std::string solver, const Rigid_Mesh & mesh, const BoundaryConditions & bc,
                const Func & func, const DataPtr_Type & data):
        Problem<QRMatrix, QRFracture, MatrixType>(solver, mesh, bc, func, data) {};

    //! Get the stiffness matrix
    /*!
     * @return the stiffness matrix
     */
    virtual Matrix_Type & getStiffnessMatrix() { return this->M_A; }

    //! Assemble matrix method
    /*!
     * Build the stiffness matrix and the BCs.
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
    virtual void solve() { this->M_solver->solve();
        Eigen::saveMarket( this->M_A, "./mMatrix/A.m" );
        Eigen::saveMarket( this->M_b, "./mMatrix/b.m" );
        Eigen::saveMarket( this->M_solver->getSolution(), "./mMatrix/sol.m" ); }

    //! Destructor
    virtual ~DarcySteady() = default;
};

template <class QRMatrix, class QRFracture, typename MatrixType>
void DarcySteady<QRMatrix, QRFracture, MatrixType >::
assembleMatrix()
{
    this->M_quadrature.reset( new Quadrature(this->M_mesh, QRMatrix(), QRFracture()) );

    if(this->M_numet == Data::NumericalMethodType::FV)
    {
		StiffMatrixFV S(this->M_mesh, this->M_mesh.getCellsVector().size() 
			+ this->M_mesh.getFractureFacetsIdsVector().size(), this->M_bc);
        S.assemble();
        S.closeMatrix();
        S.showMe();
        
        this->M_A = S.getMatrix();    // Move semantic!!!
		this->M_b = S.getBCVector();  // Move semantic!!!
    }
    else if(this->M_numet == Data::NumericalMethodType::MFD)
    {
        StiffMatrixMFD S(this->M_mesh, this->M_mesh.getFacetsVector().size() + this->M_mesh.getFractureFacetsIdsVector().size()
			+ this->M_mesh.getCellsVector().size() + this->M_mesh.getFractureFacetsIdsVector().size(), this->M_bc);
		S.assemble();
		S.CompressMatrix();
		S.showMe();
		
		this->M_A = S.getMatrix();    // Move semantic!!!
		this->M_b = S.getBCVector();  // Move semantic!!!
    }
} // DarcySteady::assembleMatrix

template <class QRMatrix, class QRFracture, typename MatrixType>
void DarcySteady<QRMatrix, QRFracture, MatrixType >::
assembleVector()
{
    this->M_quadrature.reset( new Quadrature(this->M_mesh, QRMatrix(), QRFracture()) );
		
		UInt numCellsTot = this->M_mesh.getCellsVector().size() + this->M_mesh.getFractureFacetsIdsVector().size();
		Vector f( Vector::Constant( numCellsTot, 0.) );   
		
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
        f += this->M_quadrature->cellIntegrateFractures(this->M_func);
    } // if

    if ( this->M_b.size() == 0 )
    {
        this->M_b = f;
    } // if
    if(this->M_numet == Data::NumericalMethodType::FV)
    {
		this->M_b += f;
		}
    if(this->M_numet == Data::NumericalMethodType::MFD)
    {
		UInt numFacetsTot = this->M_mesh.getFacetsVector().size() + this->M_mesh.getFractureFacetsIdsVector().size();
        this->M_b.segment(numFacetsTot,numCellsTot) -= f;
    } 
} // DarcySteady::assembleVector

} // namespace FVCode3D
#endif /* DARCYSTEADY_HPP_ */
