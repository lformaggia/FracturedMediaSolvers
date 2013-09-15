 /*!
 *	@file DarcyTypeDefinitions.hpp
 *	@brief Typedefinitions for the namespace Darcy.
 *
 *	@author Francesco Della Porta 
 *
 */ 


#ifndef __DARCYTYPEDEFNITIONS_HPP__
#define __DARCYTYPEDEFNITIONS_HPP__

#include <Eigen/Sparse>

namespace Darcy
{

//! Type for real numbers
typedef double D_Real;
//! Type for integer numbers
typedef long unsigned int D_UInt;
//! Type for eigen column-major sparse matrix of D_Real
typedef Eigen::SparseMatrix<D_Real> SpMat;
//! Type for eigen vectors
typedef Eigen::VectorXd Vector;
//! Type for eigen triplets
typedef Eigen::Triplet<D_Real> Triplet;

}

#endif
