/*!
 *	@file MatrixHandler.hpp
 *	@brief Base Class for building a matrix on Rigid_Mesh by discretizing a PDE.
 *
 *	@author Francesco Della Porta 
 *
 */ 

#ifndef __MATRIXHANDLER_HPP__
#define __MATRIXHANDLER_HPP__

#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
#include "DarcyTypeDefinitions.hpp"
	
namespace Darcy{

/*!
	@enum DiscretizationType
	It is possible to choose if the Matrix has the size.
*/
enum DiscretizationType {D_Cell, D_Nodes};

/*!
	@class MatrixHandler
	This class is a base class. The aim of this class is to be a base from which derive classes that build matrix from discretizing PDE's on a Rigid_Mesh
*/
template <class T> 
class MatrixHandler {

	typedef typename T::Generic_Point Generic_Point;
	typedef typename T::Generic_Vector Generic_Vector;
	typedef DiscretizationType DType;

public:
	//! @name Constructor & Destructor
	//@{

	//! Constructor for a MatrixHandler, given a Geometry::Rigid_Mesh.
	/*!
		@param rigid_mesh A Geometry::Rigid_Mesh on which the matrix is constructed
		@param dtype Is the policy adopted in discretization: per Cell or per Nodes (default cell)
	*/
	MatrixHandler(const Geometry::Rigid_Mesh<T> &rigid_mesh, DType dtype=D_Cell);
	//! Copy-Constructor deleted
	MatrixHandler(const MatrixHandler&) = delete;
	//! Empty-Constructor deleted
	MatrixHandler() = delete;
	//! Default destructor
	virtual ~MatrixHandler() {};
	//@}

	//! @name Get Methods
	//@{
		
	//! Get Matrix (const)
	/*!
	 * @return A reference to the matrix
	 */
	SpMat& getMatrix() const
		{return *M_Matrix;}

	//! Get size (const)
	/*!
	 * @return The size N of the matrix which is a square-Matrix (N x N)
	 */
	UInt getSize() const
		{return M_size;}
	//@}

	//! @name Methods
	//@{
		
	//! assemble
	/*!
	 * @return Constructs the matrix
	 */
	virtual void assemble()=0;
	//@}

protected:
	//! A reference to a Geometry::Rigid_Mesh
	const Geometry::Rigid_Mesh<T> & M_mesh;
	//! It explains if the matrix has the size of the number of Cells or the number of Nodes
	DType M_policy;
	//! The size N of the NxN matrix
	D_UInt M_size;
	//! A unique pointer to the assembled matrix
	std::unique_ptr<SpMat> M_Matrix;
};


// --------------------   Class MatrixHandler   --------------------
 
// ==================================================
// Constructors & Destructor
// ==================================================

template<class T>
MatrixHandler<T>::MatrixHandler (const Geometry::Rigid_Mesh<T> &rigid_mesh, DType dtype):M_mesh (rigid_mesh), M_policy(dtype), 
		 M_size ((1-dtype)*(rigid_mesh.getCellsVector().size()+ rigid_mesh.getFractureFacetsIdsVector().size())+dtype*(rigid_mesh.getNodesVector().size())), M_Matrix(new SpMat(this->M_size, this->M_size))
{}

}

#endif
