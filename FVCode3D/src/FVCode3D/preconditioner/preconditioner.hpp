/*!
 * @file preconditioner.hpp
 * @brief These classes implement simple classes that perform the "inversion" of the preconditioner.
 */

// Add this macro to enable the lumping of Mblock in place of taking the diagonal part
// in BlockTr and ILU preconditioners
//#define LUMP

#ifndef __PRECONDITIONER_HPP__
#define __PRECONDITIONER_HPP__

#include <FVCode3D/core/BasicType.hpp>
#include <Eigen/LU>
#include <Eigen/SparseCholesky>
#include <unsupported/Eigen/SparseExtra>

namespace FVCode3D
{

class preconditioner;
	
//! Typedef for preconPtr_Type
/*!
 * @typedef preconPtr_Type
 * This type definition permits to handle a std::shared_ptr<preconditioner> as a preconPtr_Type.
 */
typedef std::shared_ptr<preconditioner> preconPtr_Type;

//! Class for assembling a saddle point matrix.
/*!
 * @class SaddlePointMat
 * This class implements a generic saddle point matrix. It consists of 3 block matrices
 * M, B and T. It will be build up through an object of type SaddlePoint_StiffMatrix
 * that is the assembler of the numerical method and it overloads the *(Vector) operator
 * to use the blocks matrix to compute the matrix-vector product. In this way we avoid
 * to store the whole system matrix and we use the same 3 blocks (without any copy of them)
 * to make everything we need: matrix-vector product, building and inverting the preconditioner.
 * Clearly, this class is interesting only with an iterative system solving the system.
 * Constructors, assignement, destructor and so on are all defaulted.
*/
class SaddlePointMat
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Empty-Constructor
    SaddlePointMat(): 
		isSymUndef(1){}
    
    //! Construct a saddle point matrix
    /*!
     * @param SP_Stiff A reference to the saddle point stiffness matrix
     */
	SaddlePointMat(const SpMat & Mmat,const SpMat & Bmat,const SpMat & Tmat):
		isSymUndef(1), M(Mmat), B(Bmat), T(Tmat) {}
		
	//! Construct a saddle point matrix
    /*!
     * @param Mdim M block dimension
     * @param Brow B block row
     * @param Bcol B block col
     */
	SaddlePointMat(const UInt Mdim,const UInt Brow):
		isSymUndef(1), M(Mdim,Mdim), B(Brow,Mdim), T(Brow,Brow) {}
		
	//! Copy-Constructor
    SaddlePointMat(const SaddlePointMat &) = default;
	//! Destructor
    ~SaddlePointMat() = default;
	//@}
		
    //! @name Methods
    //@{
	//! Set the saddle point matrix
    /*!
     * @param SP_Stiff A reference to the saddle point stiffness matrix
     */
    void Set(const SpMat & Mmat,const SpMat & Bmat,const SpMat & Tmat)
	{
		M = Mmat;
		B = Bmat;
		T = Tmat;
	}
	
	//! Set the flag isSymUndef
    /*!
     * @param coeff The flag to be set
     */
    void Set_isSymUndef(const int coeff)
	{
		isSymUndef = coeff;
	}
	
	//! Compress the block matrices
    /*!
     * Compress the block matrices M, B and T
     */
    void makeCompressed() 
        {
			M.makeCompressed();
			B.makeCompressed();
			T.makeCompressed();
		};
    //@}
    
    //! @name Get Methods
    //@{
    //! Get M block (read only)
    /*!
     * @return A const reference to the M block
     */
    const SpMat & getM() const
        {return M;}
        
    //! Get M block 
    /*!
     * @return A reference to the M block
     */
    SpMat & getM() 
        {return M;}
    
    //! Get B block (read only)
    /*!
     * @return A const reference to the B block
     */
    const SpMat & getB() const
        {return B;}
        
    //! Get B block 
    /*!
     * @return A reference to the B block
     */
    SpMat & getB()
        {return B;}
        
    //! Get T block (read only)
    /*!
     * @return A const reference to the T block
     */
    const SpMat & getT() const
        {return T;}
        
    //! Get T block
    /*!
     * @return A reference to the T block
     */
    SpMat & getT()
        {return T;}
    //@}

    //! @name Methods
    //@{
	//! Resize the system
    /*!
     * @param SP_Stiff A reference to the saddle point stiffness matrix
     */
    void resize(const UInt Mdim,const UInt Brow)
	{
		M.resize(Mdim,Mdim);
		B.resize(Brow,Mdim);
		T.resize(Brow,Brow);
	}
	
	//! Get the number of non zero
    /*!
     * @return the number of non zero
     */
    UInt nonZeros() { return M.nonZeros()+2*B.nonZeros()+T.nonZeros(); }
    //@}

    //! @name Operators
    //@{
	//! Overload of matrix-vector product operator using the blocks
    /*!
     * @param x A const reference to the vector
     * @return The vector resulting from the matrix-vector product 
     */
    Vector operator * (const Vector & x) const
    {
		Vector result(M.rows()+B.rows());
		result.head(M.rows()) = M*x.head(M.cols()) + B.transpose()*x.tail(B.rows());
		result.tail(B.rows()) = isSymUndef*B*x.head(M.cols()) + isSymUndef*T*x.tail(B.rows());
		return result;
	}
    //@}

private:
	//! A flag indicating if it is sym-undef or defpos-unsym 
	int       isSymUndef;
	//! The M block matrix
	SpMat      M;
	//! The B block matrix 
	SpMat      B;
	//! The T block matrix
	SpMat      T;
};

//! This utility function computes the Approximate Inverse of the Mimetic Inner Product Matrix
/*!
 * @param Mat a SaddlePointMatrix
 * @param lumping. If true I return the lumped innerproduct matrix, if not I return the diagonal part.
 * @return The Shur Complement w.r.t. w.r.t. the lower diagonal block B D^-1 B^T + T
 */
DiagMat ComputeApproximateInverseInnerProd(const SaddlePointMat & Mat, bool lumping);

//! Compute the approximate Shur complement using a diagonal approximate inverse
/*!
 * @param Mat A saddle point matrix
 * @param D a diagonal matrix that represents an approximate inverse of the inner product matrix
 * @return th approximate Shur complement
 *
 */
SpMat ComputeApproximateSchur(const SaddlePointMat & Mat, const DiagMat & D);

//! Base class for assembling a preconditioner.
/*!
 * @class preconditioner
 * This is the base class for a preconditioner. As derived classes we can have any preconditioner:
 * identity, diagonal, inexact LU factorization and so on. The goal of these classes is not the
 * assemblation of the preconditioner matrix, but solving the linear system Pz = r in the iterative
 * scheme. They act as a wrapper that allows the practical usage of preconditioner.
*/
class preconditioner
{
public:
    //! Set the preconditioner
    /*!
     * @param Mat The saddle point mat
     */
    virtual void set(const SaddlePointMat & Mat) = 0;
    
    //! @name Solve Methods
    //@{
    //! Solve the linear system Pz=r.
    /*!
     * @param r The rhs vector on what we apply the P^-1
     * Note that the return value is passed through move semantic.
     */
    virtual Vector solve(const Vector & r) const = 0;
    //@}                  	
};


//! Class for assembling an identity preconditioner
/*!
 * @class identity_preconditioner
 * This class constructs an identity preconditioner. Create and pass an object of this type
 * in the case of running an iterative solver without preconditioning.
 */
class identity_preconditioner: public preconditioner
{
public:
    //! Set the use of the lumped form for the diagonal approximation of M

    //! Set the preconditioner
    /*!
     * @param Mat The saddle point mat
     */
    void set(const SaddlePointMat &)
    {
		std::cout<<"WARNING: the system is highly ill conditioned. It is highly recommended to precondition the system."<<std::endl<<std::endl;
	}
     //! @name Solve Methods
    //@{
    //! Solve the linear system Pz=r
    /*!
     * @param r The rhs vector on what we apply the P^-1
     */
    Vector solve(const Vector & r) const
		{ return r; };
    //@}                 
};

//! Class for assembling a diagonal preconditioner
/*!
 * @class diagonal_preconditioner
 * This class constructs a diagonal preconditioner.
 */
class diagonal_preconditioner: public preconditioner
{
public:
    //! @name Constructor & Destructor
    //@{
    //! Construct a diagonal preconditioner.
    /*!
     * @param Mat The matrix of the system that we want to precondition
     */
	diagonal_preconditioner( const SaddlePointMat & Mat ):
		Mptr(& Mat.getM()), Tptr(& Mat.getT()) {}
	//! Copy-Constructor deleted
	/*!
	 * The reason is that this class is just a wrapper to an existing
	 * SaddePointMat. If you activate copy it will be just a shallow copy and
	 * this may confuse the user. Better avoid copy constructions
	 *
	 * @todo It is probably simpler and less error prone to delete the copy constructor
	 * and assignment operator at the level of the base class, thus invalidating
	 * those of the derived classes automatically.
	 */
    diagonal_preconditioner(const diagonal_preconditioner &) = delete;
	//! Empty-Constructor
    diagonal_preconditioner(): Mptr(nullptr), Tptr(nullptr) {}
	//! Destructor
    ~diagonal_preconditioner() = default;
	//@}

    //! Set the preconditioner
    /*!
     * @param Mat The saddle point mat
     */
    void set(const SaddlePointMat & Mat)
    {
		Mptr = & Mat.getM();
		Tptr = & Mat.getT();
	}

    //! @name Solve Methods
    //@{
    //! Solve the linear system Pz=r
    /*!
     * @param r The rhs vector on what we apply the P^-1
     */
    Vector solve(const Vector & r) const;
    //@}
    
private:            
	//! A const pointer to the inner product matrix
    const SpMat *       Mptr;    
	//! A const pointer to the T block matrix
    const SpMat *       Tptr;               
};

//! Base class for assembling a Block Triangular preconditioner.
/*!
 * @class BlockTriangular_preconditioner
 * Class that builds up a block triangolar preconditioner. The M block is
 * approximated through the diagonal part of M and the SC through the inverse
 * of the diagonal part of M. The SC system is solved through the CG method.
*/
class BlockTriangular_preconditioner: public preconditioner
{
public:
    //! @name Constructor & Destructor
    //@{
	//! No Copy-Constructor
    BlockTriangular_preconditioner(const BlockTriangular_preconditioner &) = delete;
	//! Empty-Constructor
    BlockTriangular_preconditioner(): Bptr(nullptr) {}
	//! Destructor
    ~BlockTriangular_preconditioner() = default;
	//@}
    /*!
      *  If we set lumping=true the approximation of the matrix M in the Shur complement will be made
      *  by the lumping strategy, i.e. by summing the rows. The default (lumping=false) uses the diagonal part of M
      */
       void setLumping(bool lumping)
       {
         this->lumped=lumping;
       }
       //! To verify if lumping is used instead of diagonal to approximate M
       bool lumping() const {return this->lumped;}
    //! @name Assemble Methods
    //@{
	//! Assemble the the inverse diag of M and the SC
    /*!
     * @param Mat The saddle point mat
     *
     */
       void set(const SaddlePointMat & SP) override
       {
         Bptr   = & SP.getB();
         Md_inv=ComputeApproximateInverseInnerProd(SP,lumped);
         // B D^{-1}B^T -T
         chol.compute(-ComputeApproximateSchur(SP, Md_inv));// store chowleski factor
       }
    //@}

    //! @name Solve Methods
    //@{
    //! Solve the linear system Pz=r
    /*!
     * @param r The rhs vector on what we apply the P^-1
     */
    Vector solve(const Vector & r) const;
    //@}
    
private:
	//! The B block matrix
    const SpMat *              Bptr;
    //! The inverse of the diagonal of M
    DiagMat                    Md_inv;
    //! Cholesky factorization
    Eigen::SimplicialLDLT<SpMat, Eigen::Upper> chol;
    bool lumped=false;
};

//! A Symmetric block diagonal preconditioner to be used with MINRES
/*!
 *  Following the indication of the book of Wather et al.
 *
 */
class BlockDiagonal_preconditioner: public preconditioner
{
public:
    //! @name Constructor & Destructor
    //@{
      //! No Copy-Constructor
    BlockDiagonal_preconditioner(const BlockDiagonal_preconditioner &) = delete;
      //! Empty-Constructor
    BlockDiagonal_preconditioner(): Bptr(nullptr) {}
      //! Destructor
    ~BlockDiagonal_preconditioner() = default;
      //@}
    /*!
      *  If we set lumping=true the approximation of the matrix M in the Shur complement will be made
      *  by the lumping strategy, i.e. by summing the rows. The default (lumping=false) uses the diagonal part of M
      */
       void setLumping(bool lumping)
       {
         this->lumped=lumping;
       }
       //! To verify if lumping is used instead of diagonal to approximate M
       bool lumping() const {return this->lumped;}
    //! @name Assemble Methods
    //@{
      //! Assemble the the inverse diag of M and the SC
    /*!
     * @param Mat The saddle point mat
     *
     */
       void set(const SaddlePointMat & SP) override
       {
         Bptr   = & SP.getB();
         Md_inv = ComputeApproximateInverseInnerProd(SP,lumped);
         chol.compute(-ComputeApproximateSchur(SP, Md_inv));// store chowleski factor
       }
    //@}

    //! @name Solve Methods
    //@{
    //! Solve the linear system Pz=r
    /*!
     * @param r The rhs vector on what we apply the P^-1
     */
    Vector solve(const Vector & r) const;
    //@}

private:
      //! The B block matrix
    const SpMat *              Bptr;
    //! The inverse of the diagonal of M
    DiagMat                    Md_inv;
    //! Cholesky factorization
    Eigen::SimplicialLDLT<SpMat, Eigen::Upper> chol;
    bool lumped=false;
};

//! Class for assembling a ILU preconditioner.
/*!
 * @class ILU_preconditioner
 * This class builds up a ILU preconditioner where the M block is approximated via its diagonal part
 * and the SC is approximated through the inverse of the diagonal part of M.
 * The SC linear system is solved through the Eigen CG.
 * @todo again it is better to have the version for the lumped
*/
class ILU_preconditioner: public preconditioner
{
public:
    //! @name Constructor & Destructor
    //@{
	//! No Copy-Constructor
    ILU_preconditioner(const ILU_preconditioner &) = delete;
	//! Empty-Constructor
    ILU_preconditioner(): Bptr(nullptr) {}
	//! Destructor
    ~ILU_preconditioner() = default;
	//@}
    /*!
       *  If we set lumping=true the approximation of the matrix M in the Shur complement will be made
       *  by the lumping strategy, i.e. by summing the rows. The default (lumping=false) uses the diagonal part of M
       */
        void setLumping(bool lumping)
        {
          this->lumped=lumping;
        }
        //! To verify if lumping is used instead of diagonal to approximate M
        bool lumping() const {return this->lumped;}

    //! @name Assemble Methods
    //@{
	//! Assemble the approximations of M and SC
    /*!
     * @param Mat The saddle point mat
     */
        void set(const SaddlePointMat & SP)
        {
          Bptr   = & SP.getB();
          Md_inv=ComputeApproximateInverseInnerProd(SP,lumped);
          chol.compute(-ComputeApproximateSchur(SP, Md_inv));// store chowleski factor
        }
        //@}

    //! @name Solve Methods
    //@{
    //! Solve the linear system Pz=r
    /*!
     * @param r The rhs vector on what we apply the P^-1
     */
    Vector solve(const Vector & r) const;
    //@}
    
private:
	//! The B block matrix
    const SpMat *              Bptr;
    //! The inverse of the diagonal of M
    DiagMat                    Md_inv;
    //! Cholesky factorization
    Eigen::SimplicialLDLT<SpMat, Eigen::Upper> chol;
    bool lumped=false;
};

//! Class for assembling a HSS preconditioner.
/*!
 * @class HSS_preconditioner
 * This class builds up a HSS preconditioner.
 * The linear systems are solved through the Eigen CG.
*/
class HSS_preconditioner: public preconditioner
{
public:
    //! @name Constructor & Destructor
    //@{
	//! No Copy-Constructor
    HSS_preconditioner(const HSS_preconditioner &) = delete;
	//! Empty-Constructor
    HSS_preconditioner():Ncell(0), Nfrac(0), Bptr(nullptr), alpha(alpha_default), MaxIt(MaxIt_default), tol(tol_default) {}
	//! Destructor
    ~HSS_preconditioner() = default;
	//@}
    
    //! @name Set Methods
    //@{
    
    //! Set the max it value for CG
    /*!
     * @param itmax Max it value for CG
     */
	void set_alpha(const UInt a)
		{alpha = a;}
    
    //! Set the max it value for CG
    /*!
     * @param itmax Max it value for CG
     */
	void set_MaxIt(const UInt itmax)
		{MaxIt = itmax;}
		
	//! Set the tolerance value for CG
    /*!
     * @param Tol tolerance value for CG
     */
	void set_tol(const UInt Tol)
		{tol = Tol;}
    //@}

    //! @name Assemble Methods
    //@{
	//! Assemble the approximations of M and SC
    /*!
     * @param SP_Stiff A reference to the saddle point stiffness matrix
     */
    void set(const SaddlePointMat & SP);
    //@}

    //! @name Solve Methods
    //@{
    //! Solve the linear system Pz=r
    /*!
     * @param r The rhs vector on what we apply the P^-1
     */
    Vector solve(const Vector & r) const;
    //@}
    
private:
    //! Number of cells
    UInt                       Ncell;
    //! Number of farcture facets
    UInt                       Nfrac;
	//! A  constant reference to the saddle point matrix
	const SpMat *     		   Bptr;
	//! The matrix Halpha
	SpMat                      Halpha;
    //! The CG for Halpha
    Eigen::ConjugateGradient<SpMat> cg;
    // Cholesky factorization for Talpha
    Eigen::SimplicialLDLT<SpMat, Eigen::Upper> cholT;
    //! Cholesky factorization for BBtalpha
    Eigen::SimplicialLDLT<SpMat, Eigen::Upper> cholBBt;
    //! The alpha coefficient of the scheme
    Real                       alpha;
    //! The max it for CG
    UInt                       MaxIt;
    //! The tolerance for CG
    Real                       tol;
    //! The alpha coeff of the scheme (default value)
    static constexpr Real      alpha_default = 1e-2;
    //! The max it for CG (default value)
    static constexpr UInt      MaxIt_default = 300;
    //! The tolerance for CG (default value)
    static constexpr Real      tol_default = 1e-2;
                      
};

}

#endif // __LOCAL_OPERATOR_HPP__
