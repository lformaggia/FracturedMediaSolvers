#!/bin/bash
#set -x
#module load gcc-glibc/5
#module load eigen/3.3.3
#module load suitesparse/4.5.1
#module load lapack
HERE=`pwd`
# where I install
TARGET_DIR=${HERE}/../fvcode3d-install
# Varius packages location, from module environment
CMAKE_DIR=/usr/
BLAS_LIB=/usr/lib/
EIGEN_DIR=/usr/local/include/eigen3
LAPACK_DIR=/usr/lib/
SUITESPARSE_DIR=/usr/lib/x86_64-linux-gnu
./install-devel.sh --skip-packages-installation --prefix=${TARGET_DIR} \
--with-cmake=${CMAKE_DIR} --with-blas=${BLAS_LIB} --with-lapack=${LAPACK_DIR} \
--with-eigen=${EIGEN_DIR} --with-suitesparse=${SUITESPARSE_DIR} $*


