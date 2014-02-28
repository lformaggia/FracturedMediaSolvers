/*!
 *	@file basic_type.hpp
 *	@brief Basic definition types
 */ 

#ifndef BASIC_TYPE_HPP_
#define BASIC_TYPE_HPP_

#include <stdint.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <memory>
#include <iomanip>
#include <utility>
#include <functional>
#include <tuple>

#include <Eigen/Sparse>

//! Type for real numbers
typedef double Real;

//! Type for integer numbers
typedef uint64_t UInt;

//! Generic integer data
typedef int64_t Int;

//! bit-flag with up to 8 different flags
typedef unsigned char Flag8bit;

//! bit-flag with up to 32 different flags
typedef uint32_t Flag32bit;

//! PI
const Real _PI_ = std::atan(1)*4;

//! Type for eigen column-major sparse matrix of Real
typedef Eigen::SparseMatrix<Real> SpMat;

//! Type for eigen vectors
typedef Eigen::VectorXd Vector;

//! Type for eigen triplets
typedef Eigen::Triplet<Real> Triplet;

#endif /* BASIC_TYPE_HPP_ */
