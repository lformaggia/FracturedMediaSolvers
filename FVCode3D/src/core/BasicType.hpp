/*!
 *  @file basic_type.hpp
 *  @brief Basic definition types
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
#include <string>
#include <map>

#include <Eigen/Sparse>

namespace FVCode3D
{

//! Type for real numbers
typedef double Real;

//! Type for integer numbers
typedef uint64_t UInt;

//! Generic integer data
typedef int64_t Int;

//! bit-flag with up to 8 different flags
typedef uint8_t Flag8bit;

//! bit-flag with up to 16 different flags
typedef uint16_t Flag16bit;

//! bit-flag with up to 32 different flags
typedef uint32_t Flag32bit;

//! bit-flag with up to 64 different flags
typedef uint64_t Flag64bit;

//! PI
const Real _PI_ = std::atan(1)*4;

//! Type for eigen column-major sparse matrix of Real
typedef Eigen::SparseMatrix<Real> SpMat;

//! Type for eigen vectors
typedef Eigen::VectorXd Vector;

//! Type for eigen triplets
typedef Eigen::Triplet<Real> Triplet;
class Point3D;
typedef std::function<Real (Point3D point)> Func;
} // namespace FVCode3D

#endif /* BASIC_TYPE_HPP_ */
