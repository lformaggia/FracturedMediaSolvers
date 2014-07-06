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
#include <FVCode3D/core/TypeDefinition.hpp>
#include <FVCode3D/mesh/RigidMesh.hpp>

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
		M_mesh (rigid_mesh), M_properties(M_mesh.getPropertiesMap()), M_policy(dtype),
		M_size ((1-dtype)*(rigid_mesh.getCellsVector().size()+rigid_mesh.getFractureFacetsIdsVector().size())+dtype*(rigid_mesh.getNodesVector().size())),
		M_offsetRow(0), M_offsetCol(0),
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

	//! Get Triplets (const)
	/*!
	 * @return A reference to the matrix
	 */
	const std::vector<Triplet> & getTriplets() const
		{return M_matrixElements;}

	//! Get size (const)
	/*!
	 * @return The size N of the matrix
	 */
	UInt getSize() const
		{return M_size + M_offsetRow;}
	//@}

	//! @name Methods
	//@{
		
	//! Assemble method
	/*!
	 * Assemble the matrix
	 */
	virtual void assemble()=0;

	//! Close the matrix
	/*!
	 * Fill the matrix with triplets, and clear the vector of triplets
	 */
	virtual void closeMatrix()
	{
		M_Matrix->setFromTriplets( std::begin( M_matrixElements ), std::end( M_matrixElements ) );
		M_matrixElements.clear();
	}

	//! Clear triplets
	/*!
	 * Clear the vector of triplets
	 */
	void clearTriplets()
		{ M_matrixElements.clear(); }

	//! Set offsets
	/*!
	 * @param row row offset
	 * @param col column offset
	 */
	virtual void setOffsets(const UInt row, const UInt col)
		{ M_offsetRow = row; M_offsetCol = col; M_Matrix.reset(new SpMat(M_size + M_offsetRow, M_size + M_offsetCol)); }

	//@}

protected:
	//! A reference to a Rigid_Mesh
	const Rigid_Mesh & M_mesh;
	//! A reference to a PropertiesMap
	const PropertiesMap & M_properties;
	//! It explains if the matrix has the size of the number of Cells or the number of Nodes
	DType M_policy;
	//! The size N of the NxN matrix
	UInt M_size;

	//! Row offset
	UInt M_offsetRow;
	//! Column offset
	UInt M_offsetCol;

	//! A unique pointer to the assembled matrix
	std::unique_ptr<SpMat> M_Matrix;
	//! Vector of triplets
	std::vector<Triplet> M_matrixElements;
};

} // namespace FVCode3D

#endif
