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
#include <FVCode3D/preconditioner/preconditioner.hpp>
#include <unsupported/Eigen/SparseExtra>
#include <Eigen/LU>
#include <vector>

namespace FVCode3D
{

//! Class that defines the steady-state Darcy problem
/*!
 * @class DarcySteady
 * This class defines the steady-state Darcy problem given a mesh, boundary conditions and source term.
 * It defines the assemble and the solve method.
 * The first and the second template parameter indicate the quadrature rule for the matrix and fracture respectively.
 */
template <class QRMatrix, class QRFracture>
class DarcySteady : public Problem<QRMatrix, QRFracture>
{
public:

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
        Problem<QRMatrix, QRFracture>(solver, mesh, bc, func, data) {};

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
    virtual void solve() { this->M_solver->solve(); }

    //! Destructor
    virtual ~DarcySteady() = default;
};

template <class QRMatrix, class QRFracture>
void DarcySteady<QRMatrix, QRFracture>::
assembleMatrix()
{
	const UInt numFracture  = this->M_mesh.getFractureFacetsIdsVector().size();
	const UInt numFacetsTot = this->M_mesh.getFacetsVector().size() + numFracture;
	const UInt numCell      = this->M_mesh.getCellsVector().size();
	
    this->M_quadrature.reset( new Quadrature(this->M_mesh, QRMatrix(), QRFracture()) );
    
/*	auto Algebry  = this->getAlgebry();
	auto & M_A    = std::get<0>(Algebry);
	auto & M_SP   = std::get<1>(Algebry);
	auto & M_b    = std::get<2>(Algebry);
*/
    if( this->M_numet == Data::NumericalMethodType::FV && this->M_solvPolicy == Data::SolverPolicy::Direct )
    {
		StiffMatrixFV S(this->M_mesh, numCell+numFracture, this->M_bc);
        S.assemble();
        S.closeMatrix();
        S.showMe();
        
        this->getMatrix() = S.getMatrix();    
		this->getRHS()    = S.getBCVector();  
    }
   else if( this->M_numet == Data::NumericalMethodType::MFD && this->M_solvPolicy == Data::SolverPolicy::Direct )
    {
        StiffMatrixMFD S(this->M_mesh, numFacetsTot+numCell+numFracture, this->M_bc);
		S.assemble();
		S.CompressMatrix();
		S.showMe();		
		
        this->getMatrix() = S.getMatrix();    
		this->getRHS()    = S.getBCVector();
    }
    else if( this->M_numet == Data::NumericalMethodType::MFD && this->M_solvPolicy == Data::SolverPolicy::Iterative )
    {
		SaddlePoint_StiffMatrix S(this->M_mesh, this->M_bc, numFacetsTot, numCell+numFracture, numFacetsTot);
		S.assemble();
		S.Compress();
		S.showMe();		
		
		this->getSaddlePointMatrix().Set(S.getM(),S.getB(),S.getT()); 
		this->getRHS() = S.getBCVector();
	}
} // DarcySteady::assembleMatrix

template <class QRMatrix, class QRFracture>
void DarcySteady<QRMatrix, QRFracture>::
assembleVector()
{
    this->M_quadrature.reset( new Quadrature(this->M_mesh, QRMatrix(), QRFracture()) );
	auto & M_b = this->getRHS();
		
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

    if ( M_b.size() == 0 )
    {
         M_b = f;
    } // if
    if(this->M_numet == Data::NumericalMethodType::FV)
    {
		 M_b += f;
		}
    if(this->M_numet == Data::NumericalMethodType::MFD)
    {
		UInt numFacetsTot = this->M_mesh.getFacetsVector().size() + this->M_mesh.getFractureFacetsIdsVector().size();
        M_b.segment(numFacetsTot,numCellsTot) -= f;
    } 
} // DarcySteady::assembleVector

} // namespace FVCode3D
#endif /* DARCYSTEADY_HPP_ */
