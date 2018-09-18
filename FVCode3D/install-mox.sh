#!/bin/bash
#set -x
module load gcc-glibc/5
module load eigen/3.3.3
module load suitesparse/4.5.1
module load lapack
HERE=`pwd`
# where I install
TARGET_DIR=${HERE}/../FVCODE3D_build
# Varius packages location, from module environment
CMAKE_DIR=/u/sw/pkgs/toolchains/gcc-glibc/5/base
BLAS_LIB=${mkBlasPrefix}/lib64/
EIGEN_DIR=$mkEigenInc
LAPACK_DIR=$mkLapackPrefix
SUITESPARSE_DIR=$mkSuitesparsePrefix
./install-devel.sh --skip-packages-installation --prefix=${TARGET_DIR} \
--with-cmake=${CMAKE_DIR} --with-blas=${BLAS_LIB} --with-lapack=${LAPACK_DIR} \
--with-eigen=${EIGEN_DIR} --with-suitesparse=${SUITESPARSE_DIR} $*


