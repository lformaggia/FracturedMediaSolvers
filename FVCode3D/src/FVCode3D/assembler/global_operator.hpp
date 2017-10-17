/*!
 * @file global_operator.hpp
 * @brief These classes implement global operators.
 */

#ifndef __GLOBAL_OPERATOR_HPP__
#define __GLOBAL_OPERATOR_HPP__

#include <vector>
#include <cmath>
#include <FVCode3D/mesh/RigidMesh.hpp>
#include <FVCode3D/assembler/local_operator.hpp>
#include <FVCode3D/boundaryCondition/BC.hpp>
#include <Eigen/LU>
#include <unsupported/Eigen/SparseExtra>

namespace FVCode3D
{
	
//! Base class for assembling a global operator of the problem.
/*!
 * @class global_Operator
 * This is the base class to build one of the operators of the problem (the operator
 * can be a mimetic one associated to the bulk problem or a fracture operator such
 * as the coupling conditions or the finite volume trasmissibility matrix).
 * It's an abstract class.
 */
class global_Operator
{
public:
    //! Define the type of discretization
    /*!
     * @enum dType
     * This enumerator is used to show the type of discretization of row/column
     * of the operators.
     */
	enum dType{dCell, dFacet, dFracture};

    //! @name Constructor & Destructor
    //@{
    //! Construct a global_Operator.
    /*!
     * @param rMesh The rigid mesh of the problem
     * @param Prow The row dimension type
     * @param Pcol The col dimension type
     * @param rowSize The number of row
     * @param colSize The number of col
     */
	global_Operator( const Rigid_Mesh & rMesh, dType Prow, dType Pcol, UInt rowSize, UInt colSize ):
		M_mesh(rMesh), row_policy(Prow), col_policy(Pcol), Nrow(rowSize), Ncol(colSize),
		M_matrix(new SpMat(Nrow, Ncol)){}
	//! No Copy-Constructor
    global_Operator(const global_Operator &) = delete;
	//! No Empty-Constructor
    global_Operator() = delete;
	//! Destructor
    virtual ~global_Operator() = default;
	//@}

	//! @name Get Methods
    //@{
    //! Get the row dimension type
    /*!
     * @return the row dimension type
     */
    dType get_dTypeRow() const 
		{ return row_policy; };
		
	//! Get the column dimension type
    /*!
     * @return the column dimension type
     */
    dType get_dTypeCol() const 
		{ return col_policy; };
		
	//! Get the number of rows
    /*!
     * @return the number rows
     */
    UInt get_Nrow() const 
		{ return Nrow; };
		
	//! Get the number of columns
    /*!
     * @return the number columns
     */
    UInt get_Ncol() const 
		{ return Ncol; };
	
    //! Get the operator matrix (read only version)
    /*!
     * @return A const reference to the operator matrix
     */
    const SpMat & getMatrix_readOnly() const
		{ return *M_matrix; };
		
	//! Get the operator matrix (writable version)
    /*!
     * @return A reference to the operator matrix
     */
    SpMat & getMatrix() 
		{ return *M_matrix; };
	//@}

    //! @name Methods
    //@{
    //! Compress method
    /*!
     * Convert the matrix in compressed format
     */
    void CompressMatrix()
		{ (*M_matrix).makeCompressed(); };
    
    //! Reserve method
    /*!
     * Reserve the proper space for the operator matrix
     */
    virtual void reserve_space() = 0;
    
    //! Assemble method
    /*!
     * Assemble the global operator matrix
     */
    virtual void assemble() = 0;
    //@}

protected:
	//! Constant reference to the rigid mesh
	const Rigid_Mesh         & M_mesh;
	//! The row policy
	dType                      row_policy;
	//! The col policy
	dType                      col_policy;
	//! The number of rows
	UInt                       Nrow;
	//! The number of columns
	UInt                       Ncol;
	//! The sparse matrix representing the operator
	std::unique_ptr<SpMat>     M_matrix;                	
};


//! Class for assembling a global inner product matrix.
/*!
 * @class global_InnerProduct
 * This is the classe to assemble the global inner product matrix. The strategy used
 * is an assemblation cell by cell building a local inner product matrix and then
 * assembling local contributions in the global matrix.
 * The class allows also to impose the Neumann and Dirichlet bcs imposition both on
 * the matrix and on the right hand side.
 */
class global_InnerProduct: public global_Operator
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Construct a global_InnerProduct operator
    /*!
     * @param rMesh The rigid mesh of the problem
     * @param Prow The row dimension type
     * @param Pcol The col dimension type
     * @param rowSize The number of row
     * @param colSize The number of col
     * @param BCmap The boundary conditions
     */
	global_InnerProduct( const Rigid_Mesh & rMesh, dType Prow, dType Pcol, UInt rowSize, UInt colSize, const BoundaryConditions & BCmap ):
		global_Operator(rMesh, Prow, Pcol, rowSize, colSize), M_bc(BCmap){}
	//! No Copy-Constructor
    global_InnerProduct(const global_InnerProduct &) = delete;
	//! No Empty-Constructor
    global_InnerProduct() = delete;
	//! Destructor
    ~global_InnerProduct() = default;
	//@}

	//! @name Methods
    //@{
    //! Reserve method
    /*!
     * Reserve the proper space for the inner product matrix
     */
    void reserve_space();

    //! Assemble method
    /*!
     * Assemble the global inner product matrix
     */
    void assemble();
    
    //! Impose BC method
    /*!
     * This method impose the Neumann and Dirichlet BC on the monolithic matrix
     * and on the rhs.
     * @param S The monolithic matrix of the system
     * @param rhs The right hand side vector of the system
     */
    void ImposeBC(SpMat & S, Vector & rhs);
    
    //! Impose BC method
    /*!
     * This method impose the Neumann and Dirichlet BC on the inner product matrix
     * and on the rhs of the system.
     * @param rhs The right hand side vector of the system
     */
    void ImposeBC(Vector & rhs);
    //@}
    
private:
	//! Constant reference to the boundary conditions
	const BoundaryConditions & 	M_bc;
               	
};

//! Class for assembling a global divergence matrix.
/*!
 * @class global_Div
 * This is the classe to assemble the global divergence matrix. The strategy used
 * is an assemblation cell by cell building a local divergence matrix and then
 * assembling local contributions in the global matrix.
 * The class assemble also the matrix Dt that is the divergence matrix transposed
 * and modified to take into account the Neumann bc.
 */
class global_Div: public global_Operator
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Construct a global_Div operator
    /*!
     * @param rMesh The rigid mesh of the problem
     * @param Prow The row dimension type
     * @param Pcol The col dimension type
     * @param rowSize The number of row
     * @param colSize The number of col
     * @param BCmap The boundary conditions
     */
	global_Div( const Rigid_Mesh & rMesh, dType Prow, dType Pcol, UInt rowSize, UInt colSize, const BoundaryConditions & BCmap ):
		global_Operator(rMesh, Prow, Pcol, rowSize, colSize), M_bc(BCmap), Dt_matrix(new SpMat(Ncol, Nrow)){}
	//! No Copy-Constructor
    global_Div(const global_Div &) = delete;
	//! No Empty-Constructor
    global_Div() = delete;
	//! Destructor
    ~global_Div() = default;
	//@}

	//! @name Methods
    //@{
    //! Get the Dt matrix (read only version)
    /*!
     * @return A const reference to the operator matrix
     */
    const SpMat & getDtMatrix_readOnly() const
		{ return *Dt_matrix; };
		
	//! Get the operator matrix (writable version)
    /*!
     * @return A reference to the operator matrix
     */
    SpMat & getDtMatrix() 
		{ return *Dt_matrix; };
	
    //! Reserve method
    /*!
     * Reserve the proper space for the div matrix
     */
    void reserve_space();

    //! Assemble method
    /*!
     * Assemble the global inner div matrix
     */
    void assemble();
    
    //! Compress method for the Dt matrix
    /*!
     * Convert the matrix Dt in compressed format
     */
    void CompressDtMatrix()
		{ (*Dt_matrix).makeCompressed(); };
    //@}

private:
	//! Constant reference to the boundary conditions
	const BoundaryConditions & 	M_bc;
	//! The sparse matrix representing modified trasnposed matrix
	std::unique_ptr<SpMat>     Dt_matrix;
    
};

//! Class for assembling a coupling conditions matrix.
/*!
 * @class CouplingConditions
 * This is the classe to assemble the coupling conditions matrix. To build it a fractures
 * loop is performed. In the loop also the inner product matrix is modified.
 */
class CouplingConditions: public global_Operator
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Construct a coupling conditions matrix
    /*!
     * @param rMesh The rigid mesh of the problem
     * @param Prow The row dimension type
     * @param Pcol The col dimension type
     * @param rowSize The number of row
     * @param colSize The number of col
     * @param IP_matrix A constant reference to the inner product matrix
     */
	CouplingConditions( const Rigid_Mesh & rMesh, dType Prow, dType Pcol, UInt rowSize, UInt colSize, SpMat & IP_matrix):
		global_Operator(rMesh, Prow, Pcol, rowSize, colSize), M(IP_matrix), xsi(Default_xsi){}
	//! No Copy-Constructor
    CouplingConditions(const CouplingConditions &) = delete;
	//! No Empty-Constructor
    CouplingConditions() = delete;
	//! Destructor
    ~CouplingConditions() = default;
	//@}

	//! @name Methods
    //@{
    //! Get transpose method
    /*!
     * @return the tranpose of the coupling conditions matrix C
     */
    SpMat getTranspose() const
		{ return (*M_matrix).transpose(); };
		
	//! Get xsi parameter
    /*!
     * @return the xsi parameter of the coupling conditions 
     */
    const Real & get_xsi() const
		{ return xsi; };
		
	//! Set xsi parameter
    /*!
     * Set the xsi parameter of the coupling conditions 
     */
    void Set_xsi(const Real & xsiToSet)
		{ xsi = xsiToSet; };
    
    //! Reserve method
    /*!
     * Reserve the proper space for the coupling conditions matrix C
     */
	void reserve_space();

    //! Assemble method
    /*!
     * Assemble the coupling conditions matrix
     */
    void assemble();
    //@}

private:
	//! A constant reference to the inner product matrix
	SpMat                        & M;
	//! The parameter defining the coupling confition, it must be less or equal 1 and greater then 0
	Real                           xsi;
	//! The default parameter defining the coupling confition, it must be less or equal 1 and greater then 0
	static constexpr Real          Default_xsi = 1.;
	
};

//! Class for assembling a coupling conditions matrix.
/*!
 * @class FluxOperator
 * This is the classe to assemble the finite volume transmissibility matrix (changed in sign)
 * in fracture. To do it a fracture loop is performed; to take into account 
 * the fractures intersection the "star-delta" transformation is employed.
 * This class allow also the imposition of bcs on fractures.
 */
class FluxOperator: public global_Operator
{
	
//! Typedef for std::pair<UInt,UInt>
/*!
* @typedef Fracture_Juncture
* This type definition permits to treat std::pair<UInt,UInt> as a Fracture_Juncture.
*/
typedef std::pair<UInt,UInt> Fracture_Juncture;
//! Typedef for Rigid_Mesh::Edge_ID
/*!
* @typedef Edge_ID
* This type definition permits to treat Rigid_Mesh::Edge_ID as a Edge_ID.
*/
typedef Rigid_Mesh::Edge_ID Edge_ID;

public:
    //! @name Constructor & Destructor
    //@{
    //! Construct a FluxOperator matrix.
    /*!
     * @param rMesh The rigid mesh of the problem
     * @param Prow The row dimension type
     * @param Pcol The col dimension type
     * @param rowSize The number of row
     * @param colSize The number of col
     */
	FluxOperator( const Rigid_Mesh & rMesh, dType Prow, dType Pcol, UInt rowSize, UInt colSize,  const BoundaryConditions & BCmap):
		global_Operator(rMesh, Prow, Pcol, rowSize, colSize), M_bc(BCmap){}
	//! No Copy-Constructor
    FluxOperator(const FluxOperator &) = delete;
	//! No Empty-Constructor
    FluxOperator() = delete;
	//! Destructor
    ~FluxOperator() = default;
	//@}

    //! @name Alpha Methods
    //@{
    //! Border center
    /*!
     * @param fj the Id of the juncture of two Fracture_Facet in 3D
     * @return The center of the juncture between two Fracure_Facet
     */
    Point3D getBorderCenter(Fracture_Juncture fj) const;
    
    //! It is called by the method assemble() and it computes the coefficient alpha
    /*!
     * @param facetId the Id of a Facet
     * @param edge A pointer to a Rigid_mesh::Edge_ID
     * @return The computed coefficient alpha
     */
    Real findAlpha (const UInt & facetId, const Edge_ID * edge) const;

    //! It is called by the method assemble() and it computes the coefficient alpha in the case of Dirichlet BC
    /*!
     * @param facetId the Id of a Facet
     * @param edge A pointer to a Rigid_mesh::Edge_ID
     * @return The computed coefficient alpha
     */
    Real findDirichletAlpha (const UInt & facetId, const Edge_ID * edge) const;

    //! It is called by the method assemble() and it computes the coefficient alpha in the case of a fracture in 3D
    /*!
     * @param fj is a Fracture_Juncture
     * @param n_Id The Id of the Fracture_Facet
     * @return The computed coefficient alpha
     */
    Real findFracturesAlpha (const Fracture_Juncture fj, const UInt n_Id) const;
    //@}
    
    //! @name Assemble Methods
    //@{
    //! Reserve method
    /*!
     * Reserve the proper space for the coupling conditions matrix C
     */
	void reserve_space();

    //! Assemble method
    /*!
     * Assemble the coupling conditions matrix
     */
    void assemble();
    
    //! Impose fractures BC method
    /*!
     * This method impose the Neumann and Dirichlet BC on fractures
     * modifying the monolithic matrix and the rhs.
     * @param S The monolithic matrix of the system
     * @param rhs The right hand side vector of the system
     */
    void ImposeBConFractures(SpMat & S, Vector & rhs);
    
    //! Impose fractures BC method
    /*!
     * This method impose the Neumann and Dirichlet BC on fractures
     * modifying the transmissibility matrix and the rhs of the system.
     * @param rhs The right hand side vector of the system
     */
    void ImposeBConFractures(Vector & rhs);
    //@}

private:
	//! Constant reference to the boundary conditions
	const BoundaryConditions & 	M_bc;
	
};


//! Class for assembling a global bulk builder.
/*!
 * @class global_BulkBuilder
 * This is the classe to assemble the global bulk builder. It's useful using a builder 
 * of this type to assemble efficiently the bulk matrices. With this builder we can assemble
 * M, B, Dt matrices jointly with only one cells loop. This builder allow both to build up
 * M, B, Dt separately or to build up these matrices in the monolithic matrix of the system
 */
class global_BulkBuilder
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Construct a global_BulkBuilder operator.
    /*!
     * @param rMesh A constant reference to the mesh
     * @param BCmap A constant reference to the boundary conditions
     * @param M_matrix A reference to the inner product matrix
     * @param B_matrix A reference to the divergence matrix
     * @param Dt_matrix A reference to the transposed modified divergence matrix
     * @param S_matrix A reference to the monolithic matrix of the system
     */
	global_BulkBuilder( const Rigid_Mesh & rMesh, const BoundaryConditions & BCmap, SpMat & M_matrix, SpMat & B_matrix,
		SpMat & Dt_matrix, SpMat & S_matrix ):
		M_mesh(rMesh), M_bc(BCmap), M(M_matrix), B(B_matrix), Dt(Dt_matrix), S(S_matrix){}
	//! No Copy-Constructor
    global_BulkBuilder(const global_BulkBuilder &) = delete;
	//! No Empty-Constructor
    global_BulkBuilder() = delete;
	//! Destructor
    ~global_BulkBuilder() = default;
	//@}

	//! @name Methods
    //@{
    //! Reserve method for the monolithic matrix.
    /*!
     * Reserve the proper space for the monolithic matrix S
     */
    void reserve_space_Monolithic();
    
    //! Reserve method for the bulk matrices.
    /*!
     * Reserve the proper space for the M, B, Dt matrices separately
     */
    void reserve_space();

    //! Assemble method for the monolithic matrix.
    /*!
     * Assemble the M, B, Dt matrices in the monolithic matrix S
     */
    void assemble_Monolithic();
    
    //! Assemble method for the bulk matrices.
    /*!
     * Assemble the M, B, Dt matrices separately
     */
    void assemble();
    //@}

private:
	//! Constant reference to the rigid mesh
	const Rigid_Mesh              & M_mesh;
	//! Constant reference to the boundary conditions
	const BoundaryConditions      &	M_bc;
	//! A reference to the inner product matrix
	SpMat                         & M;
	//! A reference to the divergence matrix
	SpMat                         & B;
	//! A reference to the transposed divergence modified matrix
	SpMat                         & Dt;
	//! A reference to the monolithic matrix of the system
	SpMat                         & S; 
};

//! Class for assembling a global bulk builder.
/*!
 * @class FractureBuilder
 * This is the classe to assemble the fracture builder. It's useful using a builder 
 * of this type to assemble efficiently the fracture matrices. With this builder we can assemble
 * T and C matrices and modifying M jointly with only one fractures loop. 
 * This builder allow both to build up T and C and modifying M separately or
 * to build up these matrices in the monolithic matrix of the system.
 */
class FractureBuilder
{
	
//! Typedef for std::pair<UInt,UInt>
/*!
* @typedef Fracture_Juncture
* This type definition permits to treat std::pair<UInt,UInt> as a Fracture_Juncture.
*/
typedef std::pair<UInt,UInt> Fracture_Juncture;

public:
    //! @name Constructor & Destructor
    //@{
    //! Construct a FractureBuilder operator.
    /*!
     * @param rMesh A constant reference to the mesh
     * @param C_matrix A reference to the coupling conditions matrix
     * @param T_matrix A reference to the transmissibility matrix
     * @param M_matrix A reference to the inner product matrix
     * @param S_matrix A reference to the monolithic matrix of the system
     */
	FractureBuilder( const Rigid_Mesh & rMesh, SpMat & C_matrix, FluxOperator & fo, SpMat & M_matrix, SpMat & S_matrix, const Real & xsiOfC ):
		M_mesh(rMesh), C(C_matrix), FO(fo), M(M_matrix), S(S_matrix), xsi(xsiOfC){}
	//! No Copy-Constructor
    FractureBuilder(const FractureBuilder &) = delete;
	//! No Empty-Constructor
    FractureBuilder() = delete;
	//! Destructor
    ~FractureBuilder() = default;
	//@}

	//! @name Methods
    //@{
    //! Reserve method for the monolithic matrix.
    /*!
     * Reserve the proper space for the monolithic matrix S
     */
    void reserve_space_Monolithic();
    
    //! Reserve method for the fracture matrices.
    /*!
     * Reserve the proper space for the C, T and M matrices
     */
    void reserve_space();

    //! Assemble method for the monolithic matrix.
    /*!
     * Assemble the C, T and M matrices in the monolithic matrix S
     */
    void assemble_Monolithic();
    
    //! Assemble method for the bulk matrices.
    /*!
     * Assemble the C, T and M matrices separately
     */
    void assemble();
    //@}

private:
	//! Constant reference to the rigid mesh
	const Rigid_Mesh              & M_mesh;
	//! A reference to the coupling conditions matrix
	SpMat                         & C;
	//! A reference to the transmissibility matrix
	FluxOperator                  & FO;
	//! A reference to the inner product matrix
	SpMat                         & M;
	//! A reference to the monolithic matrix of the system
	SpMat                         & S;
	//! The xsi parameter of coupling conditions
	const Real                    & xsi; 
};


}

#endif // __GLOBAL_OPERATOR_HPP__
