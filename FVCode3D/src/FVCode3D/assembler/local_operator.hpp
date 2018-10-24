/*!
 * @file local_operator.hpp
 * @brief These classes implement local mimetic operators.
 */

#ifndef __LOCAL_OPERATOR_HPP__
#define __LOCAL_OPERATOR_HPP__

#include <vector>
#include <cmath>
#include <FVCode3D/mesh/RigidMesh.hpp>
#include <FVCode3D/core/BasicType.hpp>
#include <Eigen/LU>
#include <unsupported/Eigen/SparseExtra>

namespace FVCode3D
{

//! Base class for assembling a local mimetic operator.
/*!
 * @class local_MimeticOperator
 * This is the base class for a local mimetic operator (such as the the 
 * local face inner product and the local mimetic divergence).
 * It's an abstract class.
 */
class local_MimeticOperator
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Construct a local_MimeticOperator, given a cell.
    /*!
     * @param rMesh The rigid mesh of the problem
     * @param cell The cell of the local_MimeticOperator
     */
	local_MimeticOperator( const Rigid_Mesh & rMesh, const Rigid_Mesh::Cell & cell ):
		pMesh(rMesh), cellp(cell){}
	//! No Copy-Constructor
    local_MimeticOperator(const local_MimeticOperator &) = delete;
	//! No Empty-Constructor
    local_MimeticOperator() = delete;
	//! Destructor
    virtual ~local_MimeticOperator() = default;
	//@}

	//! @name Get Methods
    //@{
    //! Get the number of dofs of the local operator
    /*!
     * @return the number of dofs
     */
    UInt getDofs() const 
		{ return cellp.getFacetsIds().size(); };
	
	//! Get the cell to which the local inner product is associated
    /*!
     * @return a const reference to the cell
     */
    const Rigid_Mesh::Cell & getCell() const 
		{ return cellp; };
    //@}

    //! @name Methods
    //@{
    //! Assemble method
    /*!
     * Assemble the local_MimeticOperator
     */
    virtual void assemble() = 0;
    //@}

protected:
	//! Constant reference to the rigid mesh
	const Rigid_Mesh                  & pMesh;          
	//! The cell
	const Rigid_Mesh::Cell            & cellp;                    	
};


//! Class for assembling a local inner product matrix
/*!
 * @class local_InnerProduct
 * This class constructs the local inner product matrix given a cell.
 * The adopted scheme is based on building a consistency term Mp0 and a stability term Mp1
 * and then summing up the two contribution to obtain Mp.
 * For the stability term we use the scaled-orthogonal projector.
 * The class stores all the matrices needed to build Mp (that are Np,Rp,Mp0,Mp1) because
 * it may be of interest studying them (especially the consistency/stability terms).
 * Thanks to the free_unusefulSpace() methods, if one is interestd only in assembling
 * the local inner product matrix Mp it's possible to free the space occupied
 * by the intermediate matrices Np,Rp,Mp0,Mp1. The class can be used also to assemble
 * the approximation of the inverse of matrix Mp^-1 that is an exact inverse only of
 * the consistency part Mp0^-1.
 */
class local_InnerProduct: public local_MimeticOperator
{
public:
    //! @name Constructor & Destructor
    //@{

    //! Construct a local_InnerProduct, given a cell.
    /*!
     * @param rMesh The rigid mesh of the problem
     * @param cell The cell of the localInnerProduct
     * @param localK The permability associated to the cell
     */
	local_InnerProduct( const Rigid_Mesh & rMesh, const Rigid_Mesh::Cell & cell ):
		local_MimeticOperator(rMesh, cell){}
	//! No Copy-Constructor
    local_InnerProduct(const local_InnerProduct &) = delete;
	//! No Empty-Constructor
    local_InnerProduct() = delete;
	//! Destructor
    ~local_InnerProduct() = default;
	//@}

	//! @name Get Methods
    //@{
    //! Get the Rp matrix
    /*!
     * @return a reference to the Rp matrix
     */
    const Mat3 & getRp() const 
		{ return Rp; };
	
	//! Get the Np matrix
    /*!
     * @return a reference to the Np matrix
     */
    const Mat3 & getNp() const 
		{ return Np; };
		
	//! Get the Mp0 matrix
    /*!
     * @return a reference to the Mp0 matrix
     */
    const Mat & getMp0() const 
		{ return Mp0; };
	
	//! Get the Mp1 matrix
    /*!
     * @return a reference to the Mp1 matrix
     */
    const Mat & getMp1() const 
		{ return Mp1; };
	
	//! Get the Mp matrix
    /*!
     * @return a reference to the Mp matrix
     */
    const Mat & getMp() const 
		{ return Mp; };
    //@}

    //! @name Methods
    //@{
    //! Free method
    /*!
     * Free the unuseful space after assembling Mp
     */
    void free_unusefulSpace()
		{ 
			Np.resize(0,3);
			Rp.resize(0,3);
			Mp0.resize(0,0);
			Mp1.resize(0,0);
		};
    
    //! Assemble method
    /*!
     * Assemble the local_InnerProduct
     */
    void assemble();
    
    //! Inverse assemble method
     /*!
     * Assemble the inverse local_InnerProduct
     */
    void assemble_inv();
    //@}

friend class local_builder;

private:                     
	//! Facet normals matrix, necessary to build Mp
    Mat3                     Np;            
    //! Centroid * area matrix, necessary to build Mp
    Mat3                     Rp;            
    //! Consistency component Mp0
    Mat                      Mp0;           
    //! Stability component Mp1
    Mat                      Mp1;           
    //! Local Mp matrix for internal product
    Mat                      Mp;            
//	//! Scalar factor in Mp1 expression
	static constexpr Real    gamma = 2.;
};

//! Class for assembling a local div matrix
/*!
 * @class local_Div
 * This class constructs the local div matrix given a cell.
 * The adopted scheme is very simple and computes for every face of the mesh its signed facet
 * and builds up a cell vector that the local div matrix.
 */
class local_Div: public local_MimeticOperator
{
public:
    //! @name Constructor & Destructor
    //@{

    //! Construct a local_Div, given a cell.
    /*!
     * @param rMesh The rigid mesh of the problem
     * @param cell The cell of the localDiv
     */
	local_Div( const Rigid_Mesh & rMesh, const Rigid_Mesh::Cell & cell ):
		local_MimeticOperator(rMesh, cell){}
	//! No Copy-Constructor
    local_Div(const local_Div &) = delete;
	//! No Empty-Constructor
    local_Div() = delete;
	//! Destructor
    ~local_Div() = default;
	//@}

	//! @name Get Methods
    //@{
    //! Get the Bp matrix
    /*!
     * @return a reference to the Bp matrix
     */
    const std::vector<Real> & getBp() const 
		{ return Bp; };
	//@}

    //! @name Methods
    //@{
    //! Assemble method
    /*!
     * Assemble the local_Div
     */
    void assemble();
    //@}

friend class local_builder;

private:            
	//! Vector of signed factes of the cell
	std::vector<Real> Bp;     	
};


//! Class for construct a local builder
/*!
 * @class local_builder
 * This class is useful to assemble jointly a local inner product
 * and a local divergence operator to optimize the procedure.
 */
class local_builder
{
public:
    //! @name Constructor & Destructor
    //@{

    //! Construct a local_builder, given a cell
    /*!
     * @param InnP The local inner product
     * @param DivOp The local div operator
     * @param cell The cell of the localDiv
     */
	local_builder( local_InnerProduct & InnP, local_Div & DivOp, const Rigid_Mesh::Cell & cell ):
		IP(InnP), Div(DivOp), cel(cell){}
	//! No Copy-Constructor
    local_builder(const local_builder &) = delete;
	//! No Empty-Constructor
    local_builder() = delete;
	//! Destructor
    ~local_builder() = default;
	//@}

    //! @name Methods
    //@{
    //! Assemble method
    /*!
     * Assemble the lcoal_InnerProduct and the local_Div
     */
    void build();
    //@}

private:            
	//! The local inner product that we have to build up
    local_InnerProduct          & IP;            
    //! The local div operaotor that we have to build up
    local_Div                   & Div;            
    //! The cell
    const Rigid_Mesh::Cell      & cel;               	
};
	
}

#endif // __LOCAL_OPERATOR_HPP__
