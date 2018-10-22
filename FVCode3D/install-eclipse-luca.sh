#!/bin/bash
# where I install
HERE=`pwd`
TARGET_DIR=${HERE}/../eclipse-project
# Varius packages location, from module environment
CMAKE_DIR=/usr/
BLAS_LIB=/usr/lib/
EIGEN_DIR=/usr/local/include/eigen3b
LAPACK_DIR=/usr/lib/
SUITESPARSE_DIR=/usr/lib/x86_64-linux-gnu
/bin/bash ./install-eclipse-project.sh --skip-packages-installation --prefix=${TARGET_DIR} \
--with-cmake=${CMAKE_DIR} --with-blas=${BLAS_LIB} --with-lapack=${LAPACK_DIR} \
--with-eigen=${EIGEN_DIR} --leave-build-dir --with-suitesparse=${SUITESPARSE_DIR} $*


