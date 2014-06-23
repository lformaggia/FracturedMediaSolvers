/*!
 *	@file MatrixHandler.hpp
 *	@brief Base Class for building a matrix by discretizing a PDE with a finite volume method.
 */

#ifndef __MATRIXHANDLER_HPP__
#define __MATRIXHANDLER_HPP__

#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
#include "core/TypeDefinition.hpp"
#include "mesh/RigidMesh.hpp"

namespace FVCode3D
{

//! Type for matrix size
/*!
	@enum DiscretizationType
	It is possible to choose the size of the matrix.
*/
enum DiscretizationType {D_Cell, D_Nodes};

//! Base class for assembling a matrix
/*!
	@class MatrixHandler
	This class is a base class used as a model for the derived classes (such as the stiffness matrix and mass matrix).
	It implements a square-Matrix (N x N). It needs a Rigid_Mesh to be built.
	Abstract class.
*/
class MatrixHandler
{
public:

    //! Typedef for DiscretizationType
    /*!
        @typedef DType
        This type definition permits to treat DiscretizationType as a DType.
    */
	typedef DiscretizationType DType;

public:
	//! @name Constructor & Destructor
	//@{

	//! Construct a MatrixHandler, given a Rigid_Mesh.
	/*!
		@param rigid_mesh A Rigid_Mesh used to build the matrix
		@param dtype It is the policy adopted for the discretization: per Cell or per Nodes (default Cell)
	*/
	MatrixHandler(const Rigid_Mesh & rigid_mesh, DType dtype = D_Cell):
		M_mesh (rigid_mesh), M_policy(dtype),
		M_size ((1-dtype)*(rigid_mesh.getCellsVector().size()+rigid_mesh.getFractureFacetsIdsVector().size())+dtype*(rigid_mesh.getNodesVector().size())),
		M_Matrix(new SpMat(this->M_size, this->M_size)) {}
	//! No Copy-Constructor
	MatrixHandler(const MatrixHandler &) = delete;
	//! No Empty-Constructor
	MatrixHandler() = delete;
	//! Destructor
	virtual ~MatrixHandler() {};
	//@}

	//! @name Get Methods
	//@{
		
	//! Get Matrix (const)
	/*!
	 * @return A reference to the matrix
	 */
	SpMat & getMatrix() const
		{return *M_Matrix;}

	//! Get size (const)
	/*!
	 * @return The size N of the matrix
	 */
	UInt getSize() const
		{return M_size;}
	//@}

	//! @name Methods
	//@{
		
	//! Assemble method
	/*!
	 * @return Assemble the matrix
	 */
	virtual void assemble()=0;
	//@}

protected:
	//! A reference to a Rigid_Mesh
	const Rigid_Mesh & M_mesh;
	//! It explains if the matrix has the size of the number of Cells or the number of Nodes
	DType M_policy;
	//! The size N of the NxN matrix
	UInt M_size;
	//! A unique pointer to the assembled matrix
	std::unique_ptr<SpMat> M_Matrix;
};

} // namespace FVCode3D

#endif
