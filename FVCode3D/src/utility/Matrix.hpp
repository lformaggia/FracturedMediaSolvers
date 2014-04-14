/*!
 *	@file Matrix.hpp
 *	@brief Classes to handle matrices.
 */

#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include "core/TypeDefinition.hpp"

//! Enumerator used to select how to store the matrix
/*!
 * @enum MatrixOptions
 * This enumerator permits to decide how to store a matrix.
 */
enum MatrixOptions {
  Store    			= 0x01,
  SaveOnFile		= 0x02,
};

//! Class that implements the less operator for two matrix indices
/*!
 * @class lessOnMatrixIndices
 * Functor that, given two pairs of indices, e.g. (i,j) and (k,l),
 * returns if the first pair is less than the second one by using the lexicographic order.
 */
class lessOnMatrixIndices{
public:
	bool operator()(const std::pair<UInt, UInt> & a, const std::pair<UInt, UInt> & b) const{
		if (a.first < b.second)
			return true;
		if (a.first == b.first)
			if(a.second < b.second)
				return true;
		return false;
	};
};

//! Class that handles a sparse matrix.
/*!
 * @class SparseMatrix
 * This class allows to build a sparse matrix in COO format (row order).
 */
class SparseMatrix{
public:

	//! Default constructor
	/*!
	 * As default, it does not save the file.
	 */
	SparseMatrix();

	//! Constructor
	/*!
	 * As default, it saves the file.
	 * @param options allows to select how to store the matrix
	 * @param filename name of the file to save
	 */
	SparseMatrix(const Flag8bit options, const std::string filename);

	//! Get the filename
	/*!
	 * @return a string that contains the filename
	 */
	std::string getFilename() const
			{ return M_filename; };

	//! Get the value at (i,j)
	/*!
	 * @param i row
	 * @param j column
	 * @return the value associated with (i,j)
	 */
	Real getValue(const UInt i, const UInt j) const;

	//! Get the value at (i,j)
	/*!
	 * @param ij pair that contains the (row,column)
	 * @return the value associated with (row,column)
	 */
	Real getValue(const std::pair<UInt, UInt> ij) const;

	//! Number of non zeros
	/*!
	 * @returnthe number of non-zeros
	 */
	UInt nonZeros() const;

	//! Insert an element into the matrix
	/*!
	 * @param i row
	 * @param j column
	 * @param value value at (i,j)
	 */
	void fill(const UInt i, const UInt j, const Real value);

	//! Insert an element into the matrix
	/*!
	 * @param ij pair (row,column)
	 * @param value value at (row,column)
	 */
	void fill(const std::pair<UInt, UInt> ij, const Real value);

	//! Destructor
	~SparseMatrix();

private:

	//! Matrix
	std::map<std::pair<UInt, UInt>, Real, lessOnMatrixIndices > M_matrix;
	//! Number of non-zeros
	UInt M_nonZeros;

	//! Selector of the options
	Flag8bit M_options;

	//! Filename
	std::string M_filename;
	//! Output stream
	std::ofstream * M_file;

};

#endif /* MATRIX_HPP_ */
