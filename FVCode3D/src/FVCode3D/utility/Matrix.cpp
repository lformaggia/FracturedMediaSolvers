/*
 * Matrix.cpp
 *
 *  Created on: 01/feb/2014
 *      Author: viskius
 */

#include <FVCode3D/utility/Matrix.hpp>

namespace FVCode3D
{

SparseMatrix::SparseMatrix():
	M_nonZeros(0), M_options(Store), M_file(0) {}

SparseMatrix::SparseMatrix(const Flag8bit options, const std::string filename):
	M_nonZeros(0), M_options(options | SaveOnFile), M_filename(filename), M_file(0)
{
	M_file = new std::ofstream;
	M_file->open(M_filename.c_str(), std::ios_base::out);
}

Real SparseMatrix::getValue(const UInt i, const UInt j) const
{
	return getValue(std::make_pair(i,j));
}

Real SparseMatrix::getValue(const std::pair<UInt, UInt> ij) const
{
	if(M_matrix.size() == 0)
		return 0.;

	std::map<std::pair<UInt, UInt>, Real, lessOnMatrixIndices>::const_iterator it = M_matrix.find(ij);

	if ( it != M_matrix.end())
		return it->second;

	return 0.;
}

UInt SparseMatrix::nonZeros() const
{
	if (M_options & Store)
		return M_matrix.size();

	return M_nonZeros;
}

void SparseMatrix::fill(const UInt i, const UInt j, const Real value)
{
	fill(std::make_pair(i,j), value);
}

void SparseMatrix::fill(const std::pair<UInt, UInt> ij, const Real value)
{
	if(M_options & Store)
	{
		const std::map<std::pair<UInt, UInt>, Real, lessOnMatrixIndices>::iterator it = M_matrix.find(ij);

		if (it != M_matrix.end())
			it->second = value;
		else
			M_matrix.insert(std::make_pair(ij,value));
	}
	if(M_options & SaveOnFile)
	{
		*M_file<<ij.first<<" "<<ij.second<<" "<<value<<std::endl;
		M_nonZeros++;
	}
}

SparseMatrix::~SparseMatrix()
{
	if(M_file)
	{
		M_file->close();
		delete M_file;
	}
}

} // namespace FVCode3D
