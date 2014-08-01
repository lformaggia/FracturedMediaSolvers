/*!
 * @file RangeSupport.hpp
 * @brief Range support for Eigen::data_types
 */
#ifndef RANGESUPPORT_HPP_
#define RANGESUPPORT_HPP_

#include <Eigen/Sparse>

namespace std
{

template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
typename Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols>::Scalar*
begin(Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols>& v) {
  return v.data();
}
template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
typename Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols>::Scalar*
end(Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols>& v) {
  return v.data() + v.size();
}

template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
const typename Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols>::Scalar*
begin(const Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols>& v) {
  return v.data();
}
template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
const typename Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols>::Scalar*
end(const Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols>& v) {
  return v.data() + v.size();
}

} // namespace std

#endif /* RANGESUPPORT_HPP_ */
