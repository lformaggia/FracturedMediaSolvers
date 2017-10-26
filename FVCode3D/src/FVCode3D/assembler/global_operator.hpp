/*!
 * @file global_operator.hpp
 * @brief These classes implement global operators.
 */

#ifndef __GLOBAL_OPERATOR_HPP__
#define __GLOBAL_OPERATOR_HPP__

#include <vector>
#include <cmath>
#include <ostream>
#include <exception>
#include <FVCode3D/mesh/RigidMesh.hpp>
#include <FVCode3D/assembler/local_operator.hpp>
#include <FVCode3D/boundaryCondition/BC.hpp>
#include <Eigen/LU>
#include <unsupported/Eigen/SparseExtra>

namespace FVCode3D
{
using Eigen::Dynamic;
	
//! Define the type of discretization
/*!
* @enum dType
* This enumerator is used to show the type of discretization of row/column
* of the operators.
*/
enum class dType{dFacet, dCell, dFracture};

//! Show the permeability tensor
/*!
* @param o output stream
* @param d the dType that I want to print out
* @return output stream
*/
inline std::ostream & operator<<(std::ostream & o, const dType & d)
{
	switch(d)
	{
	case dType::dFacet    :
		o<<"Facet type";
		break;
	case dType::dCell     :
		o<<"Cell type";
		break;
	case dType::dFracture :
		o<<"Fracture type";	
	}		
	return o;
};
		
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
    //! Convert the matrix in compressed format
    /*!
     * Convert the matrix in compressed format using the proper Eigen functionality
     */
    void CompressMatrix()
		{ (*M_matrix).makeCompressed(); };
		
	//! Show basic matrix information
    /*!
     * Show basic matrix information
     */
	virtual void ShowMe() const
	{ 
		std::cout<<"Discretizzazione righe : "<<row_policy<<std::endl;
		std::cout<<"Numero righe : "<<Nrow<<std::endl;
		std::cout<<"Discretizzazione colonne : "<<col_policy<<std::endl;
		std::cout<<"Numero colonne : "<<Ncol<<std::endl<<std::endl;
	};
	//@}
  
	//! @name Assemble Methods
    //@{
    //! Reserve space for the global operator matrix.
    /*!
     * Reserve the proper space for the operator matrix computing it and then using the Eigen reserve() method.
     * Pure virtual method.
     */
    virtual void reserve_space() = 0;
    
    //! Assemble the global operator matrix
    /*!
     * Assemble the global operator matrix using Eigen insert/coeffRef methods.
     * Pure virtual method.
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
	//! Show matrix method
    /*!
     * Show basic matrix information
     */
    void ShowMe() const
    {
		std::cout<<"Basic information of inner prouct matrix M:"<<std::endl;
		global_Operator::ShowMe();
    };
	//@}

	//! @name Assemble Methods
    //@{
    //! Reserve space for the inner product matrix.
    /*!
     * Reserve the proper space for the inner product matrix computing it and then using the Eigen reserve() method.
     */
    void reserve_space();
    
    //! Assemble the local face contributions in the global system matrix 
    /*!
     * @param iloc The local facet Id
     * @param Mp The local inner product matrix
     * @param cell The actuak cell
     * @param S The system matrix 
     */
    void assembleFace(const UInt & iloc, const Eigen::Matrix<Real,Dynamic,Dynamic> & Mp,
		const Rigid_Mesh::Cell & cell, SpMat & S);
		
	//! Assemble the local face contributions in the inner product matrix 
    /*!
     * @param iloc The local facet Id
     * @param Mp The local inner product matrix
     * @param cell The actual cell
     */
    void assembleFace(const UInt & iloc, const Eigen::Matrix<Real,Dynamic,Dynamic> & Mp,
		const Rigid_Mesh::Cell & cell);

    //! Assemble the innner product matrix
    /*!
     * Assemble the inner product matrix using Eigen insert/coeffRef methods.
     */
    void assemble();
    
    //! Impose BCs on the monolithic matrix S and on the rhs
    /*!
     * @param S The monolithic matrix of the system
     * @param rhs The right hand side vector of the system
     */
    void ImposeBC(SpMat & S, Vector & rhs);
    
    //! Impose BC method on inner product matrix and on the rhs
    /*!
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
 * The class assembles also the matrix Dt that is the divergence matrix transposed
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

	//! @name Get Methods
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
	//@}
	
	//! @name Methods
    //@{
	//! Compress method for the Dt matrix
    /*!
     * Convert the matrix Dt in compressed format
     */
    void CompressDtMatrix()
		{ (*Dt_matrix).makeCompressed(); };
	
	//! Show matrix method
    /*!
     * Show basic matrix information
     */
	void ShowMe() const
    {
		std::cout<<"Basic information of divergence matrix B:"<<std::endl;
		global_Operator::ShowMe();
    };
    //@}

	//! @name Assemble Methods
    //@{
    //! Reserve space for the divergence matrix.
    /*!
     * Reserve the proper space for the divergence matrix computing it and then using the Eigen reserve() method.
     */
    void reserve_space();
    
    //! Assemble the local face contributions in the global system matrix
    /*!
     * @param iloc The local facet Id
     * @param Mp The local inner product matrix
     * @param cell The actuak cell
     * @param S The system matrix  
     */
    void assembleFace(const UInt & iloc, const std::vector<Real> & Bp,
		const Rigid_Mesh::Cell & cell, SpMat & S);
		
	//! Assemble the local face contributions in the divergence matrix
    /*!
     * @param iloc The local facet Id
     * @param Mp The local inner product matrix
     * @param cell The actuak cell 
     */
    void assembleFace(const UInt & iloc, const std::vector<Real> & Bp,
		const Rigid_Mesh::Cell & cell);

    //! Assemble the divergence matrix
    /*!
     * Assemble the divergence matrix using Eigen insert/coeffRef methods.
     */
    void assemble();
    //@}

private:
	//! Constant reference to the boundary conditions
	const BoundaryConditions & 	M_bc;
	//! The sparse matrix representing the modified trasnposed matrix
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

	//! @name Get Methods
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
	//@}
	
	//! @name Methods
	//@{
	//! Set xsi parameter
    /*!
     * Set the xsi parameter of the coupling conditions, it must be set between 0 and 1.
     * To avoid instabilities problem better 0.5 <= xsi <= 1.
     */
    void Set_xsi(const Real & xsiToSet) throw();
	
	//! Show matrix method
    /*!
     * Show basic matrix information
     */
	void ShowMe() const
    {
		std::cout<<"Basic information of coupling conditions matrix C:"<<std::endl;
		global_Operator::ShowMe();
    };
	//@}
    
	//! @name Assemble Methods
	//@{
    //! Reserve space for the coupling conditions matrix.
    /*!
     * Reserve the proper space for the coupling conditions matrix computing it and then using the Eigen reserve() method.
     */
    void reserve_space();
	
    //! Assemble the fracture facet contributions of coupling conditions matrix C and Ct in the system matrix
    /*!
     * @param facet_it The actual fracture facet
     * @param S The system matrix
     */
    void assembleFrFace(const Rigid_Mesh::Fracture_Facet & facet_it, SpMat & S);
    
    //! Assemble the fracture facet contributions of coupling conditions in the C matrix
    /*!
     * @param facet_it The actual fracture facet
     */
    void assembleFrFace(const Rigid_Mesh::Fracture_Facet & facet_it);
    
    //! Assemble the fracture facet contributions due to coupling conditions properly modifying M in the matrix system
    /*!
     * @param facet_it The actual fracture facet
     * @param S The system matrix
     */
    void assembleFrFace_onM(const Rigid_Mesh::Fracture_Facet & facet_it, SpMat & S);	
    
    //! Assemble the fracture facet contributions due to coupling conditions properly modifying M
    /*!
     * @param facet_it The actual fracture facet
     */
    void assembleFrFace_onM(const Rigid_Mesh::Fracture_Facet & facet_it);

    //! Assemble the coupling conditions matrix
    /*!
     * Assemble the coupling conditions matrix using Eigen insert/coeffRef methods.
     */
    void assemble();
    //@}

private:
	//! A constant reference to the inner product matrix
	SpMat                        & M;
	//! The parameter defining the coupling confition, it must be less or equal 1 and greater then 0
	Real                           xsi;
	//! The default parameter defining the coupling confition, it must be less or equal 1 and greater then 0
	static constexpr Real          Default_xsi = 1;
};

//! Class for assembling a trasmissibility fracture matrix.
/*!
 * @class FluxOperator
 * This is the classe to assemble the finite volume transmissibility matrix (changed in sign).
 * To do it a fracture loop is performed; to take into account 
 * the fractures intersection, the "star-delta" transformation is employed.
 * This class allow also the imposition of bcs on fractures.
 */
class FluxOperator: public global_Operator
{

public:

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
	Point3D getBorderCenter(Fracture_Juncture fj) const
	{
		return ( M_mesh.getNodesVector()[fj.first] + (M_mesh.getNodesVector()[fj.second] -
                 M_mesh.getNodesVector()[fj.first])/2. );
	};
	
    //! It is called by the method assemble() and it computes the coefficient alpha
    /*!
     * @param facetId the Id of a Facet
     * @param edge A pointer to a Rigid_mesh::Edge_ID
     * @return The computed coefficient alpha
     */
    Real findAlpha (const UInt & facetId, const Edge_ID * edge) const;

    //! It is called by the bcs method and it computes the coefficient alpha in the case of Dirichlet BC
    /*!
     * @param facetId the Id of a Facet
     * @param edge A pointer to a Rigid_mesh::Edge_ID
     * @return The computed coefficient alpha
     */
    Real findDirichletAlpha (const UInt & facetId, const Edge_ID * edge) const;

    //! It is called by the bcs method and it computes the coefficient alpha in the case of a fracture in 3D
    /*!
     * @param fj is a Fracture_Juncture
     * @param n_Id The Id of the Fracture_Facet
     * @return The computed coefficient alpha
     */
    Real findFracturesAlpha (const Fracture_Juncture fj, const UInt n_Id) const;
    //@}
    
    //! @name Methods
    //@{
    //! Show matrix method
    /*!
     * Show basic matrix information
     */
	void ShowMe() const
    {
		std::cout<<"Basic information of trasmissibility matrix T:"<<std::endl;
		global_Operator::ShowMe();
    };
    //@}
    
    //! @name Assemble Methods
    //@{
    //! Reserve space for the flux operator (trasmissibility) matrix.
    /*!
     * Reserve the proper space for the flux operatore (trasmissibility) matrix computing it and then using the Eigen reserve() method.
     */
    void reserve_space();

    //! Assemble fracture facet trasmissibilities in the system matrix
    /*!
     * Assemble the fracture facet trasmissibilities in the system matrix 
     * taking into account fracture trasformations with the "star-delta"trasformation
     * @param facet_it The actual fracture facet
     * @param S The system matrix
     */
    void assembleFrFace(const Rigid_Mesh::Fracture_Facet & facet_it, SpMat & S);
    
	//! Assemble fracture facet trasmissibilities
    /*!
     * Assemble the fracture facet trasmissibilities taking into account fractures 
     * intersections with the "star-delta" trasformation
     * @param facet_it The actual fracture facet
     */
    void assembleFrFace(const Rigid_Mesh::Fracture_Facet & facet_it);

    //! Assemble the flux operator (trasmissibility) matrix
    /*!
     * Assemble the flux operator (trasmissibility) matrix using Eigen insert/coeffRef methods.
     */
    void assemble();
    
    //! Impose fractures BC method on the monolithic matrix and on the rhs
    /*!
     * This method impose the Neumann and Dirichlet BC on fractures
     * modifying the monolithic matrix and the rhs.
     * @param S The monolithic matrix of the system
     * @param rhs The right hand side vector of the system
     */
    void ImposeBConFractures(SpMat & S, Vector & rhs);
    
    //! Impose fractures BC method on the trasmissibility matrix and on the rhs
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


//! Base class for assembling a global_Builder.
/*!
 * @class global_Builder
 * This is the base classe to assemble a global builder useful to build up the matrices of the
 * problem in an efficient way. It's an abstract class.
 */
 class global_Builder
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Construct a global_Builder.
    /*!
     * @param rMesh A constant reference to the mesh
     */
	global_Builder(const Rigid_Mesh & rMesh):
		M_mesh(rMesh){}
	//! No Copy-Constructor
    global_Builder(const global_Builder &) = delete;
	//! No Empty-Constructor
    global_Builder() = delete;
	//! Destructor
    ~global_Builder() = default;
	//@}

	//! @name Assemble Methods
    //@{
    //! Reserve the proper space for the monolithic matrix S
    /*!
     * @param S The system matrix
     */
    virtual void reserve_space(SpMat & S)=0;
    
    //! Reserve the proper space for matrices separately.
    /*!
     * Reserve the proper space for matrices separately (not assebled in the system matrix S)
     */
    virtual void reserve_space()=0;

    //! Assemble method for the monolithic matrix.
    /*!
     * Assemble the matrices in the monolithic matrix S
     * @param S The system matrix
     */
    virtual void build(SpMat & S)=0;
    
    //! Assemble method for the matrices separately.
    /*!
     * Assemble the matrices separately (not in the system matrix S)
     */
    virtual void build()=0;
    //@}

protected:
	//! Constant reference to the rigid mesh
	const Rigid_Mesh              & M_mesh;
};

//! Class for assembling a global bulk builder.
/*!
 * @class global_BulkBuilder
 * This is the classe to assemble the global bulk builder. It's useful using a builder 
 * of this type to assemble efficiently the bulk matrices. With this builder we can assemble
 * M, B, Dt matrices jointly with only one cells loop. This builder allow both to build up
 * M, B, Dt separately or to build up these matrices in the monolithic matrix of the system
 */
class global_BulkBuilder: public global_Builder
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Construct a global_BulkBuilder.
    /*!
     * @param BCmap A constant reference to the boundary conditions
     * @param M_matrix A reference to the inner product matrix
     * @param B_matrix A reference to the divergence matrix
     * @param Dt_matrix A reference to the transposed modified divergence matrix
     * @param S_matrix A reference to the monolithic matrix of the system
     */
	global_BulkBuilder(const Rigid_Mesh & rMesh, const BoundaryConditions & BCmap, global_InnerProduct & ip, global_Div & div):
		global_Builder(rMesh), M_bc(BCmap), IP(ip), Div(div){}
	//! No Copy-Constructor
    global_BulkBuilder(const global_BulkBuilder &) = delete;
	//! No Empty-Constructor
    global_BulkBuilder() = delete;
	//! Destructor
    ~global_BulkBuilder() = default;
	//@}

	//! @name Assemble Methods
    //@{
    //! Reserve space for the bulk matrices in the system matrix
    /*!
     * Reserve the proper space for M, B, Dt in the monolithic matrix S
     * @param S The system matrix
     */
    void reserve_space(SpMat & S);
    
    //! Reserve method for the bulk matrices.
    /*!
     * Reserve the proper space for the M, B, Dt matrices separately
     */
    void reserve_space();

    //! Assemble the bulk matrices in the monolithic matrix.
    /*!
     * Assemble the M, B, Dt matrices in the monolithic matrix S
     * @param S The system matrix
     */
    void build(SpMat & S);
    
    //! Assemble method for the bulk matrices.
    /*!
     * Assemble the M, B, Dt matrices separately
     */
    void build();
    //@}

private:
	//! Constant reference to the boundary conditions
	const BoundaryConditions      &	M_bc;
	//! A reference to the inner product
	global_InnerProduct           & IP;
	//! A reference to the divergence operator
	global_Div                    & Div;
};

//! Class for assembling a fracture builder.
/*!
 * @class FractureBuilder
 * This is the classe to assemble the fracture builder. It's useful using a builder 
 * of this type to assemble efficiently the fracture matrices. With this builder we can assemble
 * T and C matrices and modifying M jointly with only one fractures loop. 
 * This builder allow both to build up T and C and modifying M separately or
 * to build up these matrices in the monolithic matrix of the system.
 */
class FractureBuilder: public global_Builder
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Construct a FractureBuilder.
    /*!
     * @param C_matrix A reference to the coupling conditions matrix
     * @param T_matrix A reference to the transmissibility matrix
     * @param M_matrix A reference to the inner product matrix
     * @param S_matrix A reference to the monolithic matrix of the system
     */
	FractureBuilder( const Rigid_Mesh & rMesh, CouplingConditions & cc, FluxOperator & fo, global_InnerProduct & ip):
		global_Builder(rMesh), coupling(cc), FluxOp(fo), IP(ip){}
	//! No Copy-Constructor
    FractureBuilder(const FractureBuilder &) = delete;
	//! No Empty-Constructor
    FractureBuilder() = delete;
	//! Destructor
    ~FractureBuilder() = default;
	//@}

	//! @name Assemble Methods
    //@{
    //! Reserve space for the fractures matrices in the monolithic matrix.
    /*!
     * Reserve the proper space for C, Ct, modifying M and T in the monolithic matrix S
     * @param S The monolithic matrix
     */
    void reserve_space(SpMat & S);
    
    //! Reserve method for the fracture matrices separately.
    /*!
     * Reserve the proper space for C, T and modifying M separately
     */
    void reserve_space();

    //! Assemble the fractrues matrice in the monolithic matrix.
    /*!
     * Assemble the C, T and modify M in the monolithic matrix S
     * @param The monolithic matrix
     */
    void build(SpMat & S);
    
    //! Assemble method for the fracture matrices.
    /*!
     * Assemble C, T and modify M separately
     */
    void build();
    //@}

private:
	//! A reference to the coupling conditions matrix
	CouplingConditions            & coupling;
	//! A reference to the transmissibility matrix
	FluxOperator                  & FluxOp;
	//! A reference to the inner product matrix
	global_InnerProduct           & IP;
};

} //FVCode3D

#endif // __GLOBAL_OPERATOR_HPP__
