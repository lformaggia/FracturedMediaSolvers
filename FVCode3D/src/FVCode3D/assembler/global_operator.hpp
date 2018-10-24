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
#include <FVCode3D/core/BasicType.hpp>
#include <Eigen/LU>
#include <unsupported/Eigen/SparseExtra>

namespace FVCode3D
{
	
//! Define the type of discretization
/*!
* @enum dType
* This enumerator is used to show the type of discretization of row/column
* of the operators.
*/
enum dType{dFacet, dCell, dFracture};

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
}
		
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
	global_Operator( const Rigid_Mesh & rMesh, dType Prow, dType Pcol):
		M_mesh(rMesh), row_policy(Prow), col_policy(Pcol),
		Nrow( (row_policy==dFacet)*(M_mesh.getFacetsVector().size()+M_mesh.getFractureFacetsIdsVector().size()) +
		 (row_policy==dCell)*(M_mesh.getCellsVector().size()) + (row_policy==dFracture)*(M_mesh.getFractureFacetsIdsVector().size()) ),
		Ncol( (col_policy==dFacet)*(M_mesh.getFacetsVector().size()+M_mesh.getFractureFacetsIdsVector().size()) +
		 (col_policy==dCell)*(M_mesh.getCellsVector().size()) + (col_policy==dFracture)*(M_mesh.getFractureFacetsIdsVector().size()) ),
		M_matrix(Nrow, Ncol){}
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
    const SpMat & getMatrix() const
		{ return M_matrix; };
		
	//! Get the operator matrix (writable version)
    /*!
     * @return A reference to the operator matrix
     */
    SpMat & getMatrix() 
		{ return M_matrix; };
	//@}

    //! @name Methods
    //@{
    //! Convert the matrix in compressed format
    /*!
     * Convert the matrix in compressed format using the proper Eigen functionality
     */
    void CompressMatrix()
		{ M_matrix.makeCompressed(); };
		
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
	SpMat                      M_matrix;                	
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
	global_InnerProduct( const Rigid_Mesh & rMesh, dType Prow, dType Pcol):
		global_Operator(rMesh, Prow, Pcol) {}
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
    
    //! Assemble the local face contributions in a matrix S using the offsets
    /*!
     * @param iloc The local facet Id
     * @param Mp The local inner product matrix
     * @param cell The actual cell
     * @param S The matrix
     * @param rowoff The row offset 
     * @param coloff The col offset 
     */
    void assembleFace(const UInt iloc, const Mat & Mp, const Rigid_Mesh::Cell & cell,
		SpMat & S, const UInt rowoff, const UInt coloff) const;
		
	//! Assemble the local face contributions in the inner product matrix 
    /*!
     * @param iloc The local facet Id
     * @param Mp The local inner product matrix
     * @param cell The actual cell
     */
    void assembleFace(const UInt iloc, const Mat & Mp, const Rigid_Mesh::Cell & cell);

    //! Assemble the innner product matrix
    /*!
     * Assemble the inner product matrix using Eigen insert/coeffRef methods.
     */
    void assemble();
    //@}
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
	global_Div( const Rigid_Mesh & rMesh, dType Prow, dType Pcol):
		global_Operator(rMesh, Prow, Pcol) {}
	//! No Copy-Constructor
    global_Div(const global_Div &) = delete;
	//! No Empty-Constructor
    global_Div() = delete;
	//! Destructor
    ~global_Div() = default;
	//@}

	//! @name Get Methods
    //@{
		
    //! Get transpose method
    /*!
     * @return the tranpose of the coupling conditions matrix C
     */
    SpMat getTranspose() const
		{ return M_matrix.transpose(); };
	//@}
	
	//! @name Methods
    //@{
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
    
    //! Assemble the local face contributions of B and Bt in a matrix S using the offsets
    /*!
     * @param iloc The local facet Id
     * @param Mp The local inner product matrix
     * @param cell The actual cell
     * @param S The system matrix  
     * @param rowoff The row offset 
     * @param coloff The col offset 
     */
    void assembleFace(const UInt iloc, const std::vector<Real> & Bp, const Rigid_Mesh::Cell & cell,
		SpMat & S, const UInt rowoff, const UInt coloff, const bool transpose = true) const;
		
	//! Assemble the local face contributions in the divergence matrix
    /*!
     * @param iloc The local facet Id
     * @param Mp The local inner product matrix
     * @param cell The actual cell 
     */
    void assembleFace(const UInt iloc, const std::vector<Real> & Bp, const Rigid_Mesh::Cell & cell);

    //! Assemble the divergence matrix
    /*!
     * Assemble the divergence matrix using Eigen insert/coeffRef methods.
     */
    void assemble();
    //@}    
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
	CouplingConditions( const Rigid_Mesh & rMesh, dType Prow, dType Pcol, SpMat & IP_matrix):
		global_Operator(rMesh, Prow, Pcol), M(IP_matrix), xsi(Default_xsi){}
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
		{ return M_matrix.transpose(); };
		
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
    void Set_xsi(const Real & xsiToSet);
	
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
	
    //! Assemble the fracture facet contributions of  C and Ct in a matrix S using the offsets
    /*!
     * @param facet_it The actual fracture facet
     * @param S The system matrix
     * @param rowoff The row offset 
     * @param coloff The col offset 
     */
    void assembleFrFace(const Rigid_Mesh::Fracture_Facet & facet_it, SpMat & S, 
		const UInt rowoff, const UInt coloff, const bool transpose = true) const;
    
    //! Assemble the fracture facet contributions of coupling conditions in the C matrix
    /*!
     * @param facet_it The actual fracture facet
     */
    void assembleFrFace(const Rigid_Mesh::Fracture_Facet & facet_it);
    
    //! Assemble the fracture facet contributions due to coupling conditions properly modifying M in the matrix S using the offsets
    /*!
     * @param facet_it The actual fracture facet
     * @param S The system matrix
     * @param rowoff The row offset 
     * @param coloff The col offset 
     */
    void assembleFrFace_onM(const Rigid_Mesh::Fracture_Facet & facet_it, SpMat & S, const UInt rowoff, const UInt coloff) const;	
    
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
	static constexpr Real          Default_xsi = 0.75;
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
	FluxOperator( const Rigid_Mesh & rMesh, dType Prow, dType Pcol):
		global_Operator(rMesh, Prow, Pcol) {}
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
     * @param rowoff The row offset 
     * @param coloff The col offset 
     */
    void assembleFrFace(const Rigid_Mesh::Fracture_Facet & facet_it, SpMat & S, const UInt rowoff, const UInt coloff) const;
    
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
    //@}
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
    virtual ~global_Builder() = default;
	//@}

	//! @name Assemble Methods
    //@{
    //! Reserve the proper space for the monolithic matrix S
    /*!
     * @param S The system matrix
     */
    virtual void reserve_space(SpMat & S)=0;

    //! Assemble method for the monolithic matrix.
    /*!
     * Assemble the matrices in the monolithic matrix S
     * @param S The system matrix
     */
    virtual void build(SpMat & S)=0;
    //@}

protected:
	//! Constant reference to the rigid mesh
	const Rigid_Mesh              & M_mesh;
};

//! Class for assembling a global bulk builder.
/*!
 * @class global_BulkBuilder
 * This is the class to assemble the global bulk builder. It's useful using a builder 
 * of this type to assemble efficiently the bulk matrices. With this builder we can assemble
 * inner product and divergence matrices jointly with only one cells loop. This builder allows both to build up
 * them in the system matrix or in the saddle point blocks.
 */
class global_BulkBuilder: public global_Builder
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Construct a global_BulkBuilder.
    /*!
     * @param rMesh A constant reference to the rigid mesh
     * @param ip A reference to the inner product
     * @param div A reference to the divergence matrix
     */
	global_BulkBuilder(const Rigid_Mesh & rMesh, global_InnerProduct & ip, global_Div & div):
		global_Builder(rMesh), IP(ip), Div(div){}
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
     * Reserve the proper space for M, B, Bt in the monolithic matrix S
     * @param S The system matrix
     */
    void reserve_space(SpMat & S);
    
    //! Reserve method for the bulk matrices.
    /*!
     * Reserve the proper space for the inner product matrix and divergence in the M and B blocks
     * @param M The M block matrix
     * @param B The B block matrix
     */
    void reserve_space(SpMat & M, SpMat & B);

    //! Assemble the bulk matrices in the monolithic matrix.
    /*!
     * Assemble the M, B, Bt matrices in the monolithic matrix S
     * @param S The system matrix
     */
    void build(SpMat & S);
    
    //! Assemble method for the bulk matrices.
    /*!
     * Assemble the inner product matrix and divergence in the M and B blocks
     * @param M The M block matrix
     * @param B The B block matrix
     */
    void build(SpMat & M, SpMat & B);
    //@}

private:
	//! A reference to the inner product
	global_InnerProduct           & IP;
	//! A reference to the divergence operator
	global_Div                    & Div;
};

//! Class for assembling a fracture builder.
/*!
 * @class FractureBuilder
 * This is the class to assemble the fracture builder. It's useful using a builder 
 * of this type to assemble efficiently the fracture matrices. With this builder we can assemble
 * coupling conditions and trasmissibilities and modifying the inner product jointly with only one fractures loop. 
 * This builder allows both to build up these matrices and modifying M in the system matrix 
 * or in the saddle point block matrices.
 */
class FractureBuilder: public global_Builder
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Construct a FractureBuilder.
    /*!
     * @param rMesh A constant reference to the rigid mesh
     * @param cc A const reference to coupling conditions
     * @param ip A reference to the flux operator
     */
	FractureBuilder( const Rigid_Mesh & rMesh, CouplingConditions & cc, FluxOperator & fo):
		global_Builder(rMesh), coupling(cc), FluxOp(fo){}
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
    
    //! Reserve method for the block matrices separately.
    /*!
     * Reserve the proper space for coupling conditions and trasmissibilities in the blocks M, B and T
     * @param M The inner product matrix
     * @param B The divergence matrix
     * @param T The trasmissibility matrix
     */
    void reserve_space(SpMat & M, SpMat & B, SpMat & T);

    //! Assemble the fractrues matrice in the monolithic matrix.
    /*!
     * Assemble the C, T and modify M in the monolithic matrix S
     * @param The monolithic matrix
     */
    void build(SpMat & S);
    
    //! Assemble method for the fracture matrices.
    /*!
     * Assemble coupling conditions and trasmissibilities in the blocks M, B and T
     * @param M The inner product matrix
     * @param B The divergence matrix
     * @param T The trasmissibility matrix
     */
    void build(SpMat & M, SpMat & B, SpMat & T);
    //@}

private:
	//! A reference to the coupling conditions matrix
	CouplingConditions            & coupling;
	//! A reference to the flux operator
	FluxOperator                  & FluxOp;
};


//! Class for assembling for imposing the boundary conditions on the system.
/*!
 * @class BCimposition
 * This is the class to impose the bcs both on bulk and fractures.
 * The bulk Neumann bcs are imposed through a Nitsche approach that allows to keep the
 * symmetry of the system and (important!!!) to avoid a very unefficient prune Eigen command.
 * The bulk Dirichlet bcs are imposed directly on the rhs because they are natural bcs in our problem.
 * The fracture bcs are imposed through the typical FV approach (that is Neumann imposes the flux and
 * Dirichlet is imposed through a ghost-cell approach).
 */
class BCimposition
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Construct a BCimposition.
    /*!
     * @param rMesh A const reference to the rigid mesh
     * @param bc A const reference to the boujndary conditions
     */
	BCimposition( const Rigid_Mesh & rMesh, const BoundaryConditions & bc):
		M_mesh(rMesh), M_bc(bc), penalty(Default_penalty) {}
	//! No Copy-Constructor
    BCimposition(const BCimposition &) = delete;
	//! No Empty-Constructor
    BCimposition() = delete;
	//! Destructor
    ~BCimposition() = default;
	//@}
	
	//! @name Methods
    //@{
    //! Set penalty.
    /*!
     * Set the penalty coefficient.
     * @param gamma The penalty value to be set.
     */
    void set_penalty(const Real gamma)
		{ penalty = gamma; };
    //@}
	
	//! @name Assemble Methods
    //@{
    //! Assemble bulk bc in the system matrix.
    /*!
     * Assemble Neumann bcs through Nitshce approach and Dirichlet directly on the rhs.
     * @param S The monolithic matrix
     * @param rhs The rhs of the system
     */
    void ImposeBConBulk(SpMat & S, Vector & rhs) const; 
    
    //! Assemble bulk bc in block matrices M and B.
    /*!
     * Assemble Neumann bcs through Nitshce approach and Dirichlet directly on the rhs.
     * @param S The monolithic matrix
     * @param rhs The rhs of the system
     */
    void ImposeBConBulk(SpMat & M, SpMat & B, Vector & rhs) const; 
    
    //! Assemble bcs on fractures acting on the system matrix.
    /*!
     * Assemble the Neumann bcs imposing the flux, the Dirichlet through ghost-cells.
     * @param S The trasmissibility matrix
     * @param rhs The rhs of the system
     */
    void ImposeBConFracture(SpMat & S, Vector & rhs, FluxOperator & Fo) const;
    
    //! Assemble bcs on fractures acting on the T block matrix.
    /*!
     * Assemble the Neumann bcs imposing the flux, the Dirichlet through ghost-cells.
     * @param S The trasmissibility matrix
     * @param rhs The rhs of the system
     */
    void ImposeBConFracture_onT(SpMat & T, Vector & rhs, FluxOperator & Fo) const;
    //@}
    //@}

private:
	//! A reference to the rigid mesh
	const Rigid_Mesh            & M_mesh;
	//! A reference to the boundary conditions
	const BoundaryConditions    & M_bc;
	//! The penalty coefficient
	Real penalty;
	//! The default penalty value
	static constexpr Real Default_penalty = 1e30;
};


} //FVCode3D

#endif // __GLOBAL_OPERATOR_HPP__
