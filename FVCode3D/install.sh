#!/bin/bash
# This script installs FVCode3D and all the needed dependencies
# For a simple installation:
#   ./install.sh --prefix=<install dir>
# For a full list of option:
#   ./install.sh --help

PROJECT_NAME="FVCode3D"
LOWERCASE_PROJECT_NAME=$(echo ${PROJECT_NAME} | tr '[:upper:]' '[:lower:]')
UPPERCASE_PROJECT_NAME=$(echo ${PROJECT_NAME} | tr '[:lower:]' '[:upper:]')

PACKAGE_ARRAY=(cmake blas lapack suitesparse eigen)

###################
# Default options #
###################

ROOT_DIR=$PWD

# Set default installation directory if not already set
if [ -z "${INSTALL_DIR}" ]; then
    INSTALL_DIR=${PWD}/../install
fi

# set number of parallel builds to use
NUM_PROC=$(grep processor /proc/cpuinfo | wc -l)
if [ ${NUM_PROC} -gt 8 ]; then NUM_PROC=8; fi

# Set default build type if not already set
if [ -z "${BUILD_TYPE}" ]; then
    BUILD_TYPE="Release"
fi

INSTALL_PACKAGES=1
INSTALL_PROJECT=1
INSTALL_DOC=1
BUILD_PROJECT_TESTS="ON"
LEAVE_BUILD_DIR=""

C_COMPILER=""
CXX_COMPILER=""
FORTRAN_COMPILER=""

C_LIBRARY_DIRS=()
CXX_LIBRARY_DIRS=()
FORTRAN_LIBRARY_DIRS=()

####################
# Argument parsing #
####################
for i in "$@"
do
case $i in
    --build-dir=*)
        BUILD_DIR=$(echo $i | sed 's/[-a-zA-Z0-9]*=//')
        # Expand ~ ~+ ~- ~user and remove starting and trailing quotes
        eval BUILD_DIR=${BUILD_DIR}
        # If it's not an absolute path prefix it with ${PWD}
        [[ ${BUILD_DIR} == /* ]] || BUILD_DIR="${PWD}/${BUILD_DIR}"
    ;;
    -b=*|--build-type=*)
        ARG_BUILD_TYPE=$(echo $i | sed 's/[-a-zA-Z0-9]*=//')
        if [ "${ARG_BUILD_TYPE,,}" != "release" ] && [ "${ARG_BUILD_TYPE,,}" != "debug" ] && [ "${ARG_BUILD_TYPE,,}" != "profile" ]; then
            echo "Unrecognized build type."
            echo "Valid options are: Release, Debug and Profile"
            exit 1
        fi
        BUILD_TYPE=${ARG_BUILD_TYPE}
    ;;
    -c|--configure-only)
        CONFIGURE_ONLY=1
    ;;
    --c-compiler=*)
        C_COMPILER="$(echo $i | sed 's/[-a-zA-Z0-9]*=//')"
    ;;
    --cxx-compiler=*)
        CXX_COMPILER="$(echo $i | sed 's/[-a-zA-Z0-9]*=//')"
    ;;
    --fortran-compiler=*)
        FORTRAN_COMPILER="$(echo $i | sed 's/[-a-zA-Z0-9]*=//')"
    ;;
    --c-library-dir=*)
        C_LIBRARY_DIRS+="$(echo $i | sed 's/[-a-zA-Z0-9]*=//')"
    ;;
    --cxx-library-dir=*)
        CXX_LIBRARY_DIRS+="$(echo $i | sed 's/[-a-zA-Z0-9]*=//')"
    ;;
    --fortran-library-dir=*)
        FROTRAN_LIBRARY_DIRS+="$(echo $i | sed 's/[-a-zA-Z0-9]*=//')"
    ;;
    --install-*)
        PACKAGE_NAME=$(echo $i | sed 's/--install-//' | tr '[:upper:]' '[:lower:]')
        eval ${PACKAGE_NAME}="INSTALL" # 'INSTALL' means the package has to be installed
    ;;
    --leave-build-dir)
        LEAVE_BUILD_DIR=1
    ;;
    -p=*|--prefix=*)
        INSTALL_DIR=$(echo $i | sed 's/[-a-zA-Z0-9]*=//')
        # Expand ~ ~+ ~- ~user and remove starting and trailing quotes
        eval INSTALL_DIR=$INSTALL_DIR
        # If it's not an absolute path prefix it with ${PWD}
        [[ ${INSTALL_DIR} == /* ]] || INSTALL_DIR="${PWD}/${INSTALL_DIR}"
    ;;
    -h|--help)
        echo "Install ${PROJECT_NAME}"
        echo "Default installation directory is \"../install\" and default build type is \"Release\""
        echo "Optional arguments:"
        echo "  --build-dir=<path>          build directory;"
        echo "                              default value is \"<install_dir>/build\""
        echo "                              where <install_dir> can be specified with"
        echo "                              the --prefix option and defaults to \"../install\";"
        echo "                              unless --leave-build-dir is specified, this "
        echo "                              directory will be deleted at the end of the build"
        echo "  -b, --build-type=<type>     build ${PROJECT_NAME} in <type> mode,"
        echo "                              where <type> can be one of:"
        echo "                               - Release (default)"
        echo "                               - Debug"
        echo "                               - Profile"
        echo "  -c, --configure-only        configure ${PROJECT_NAME} and create makefiles"
        echo "                              in the build directory, but don't execute"
        echo "                              the actual build"
        echo "  --c-complier=<path>         c compiler"
        echo "  --cxx-complier=<path>       c++ compiler"
        echo "  --fortan-complier=<path>    fortran compiler"
        echo "  --c-library-dir=<path>      extra library directory to be used when compiling"
        echo "                              c programs (can be used multiple times)"
        echo "  --cxx-library-dir=<path>    extra library directory to be used when compiling"
        echo "                              c++ programs (can be used multiple times)"
        echo "  --fortran-library-dir=<path> extra library directory to be used when compiling"
        echo "                              fortran programs (can be used multiple times)"
        echo "  -h, --help                  print this help"
        echo "  --install-<package>         install <package>; <package> can be one of:"
        for package in "${PACKAGE_ARRAY[@]}"; do
            echo "                                - $package"
        done
        echo "                              (this is enabled by default for all packages)"
        echo "  --leave-build-dir           ${PROJECT_NAME} build directory won't be deleted"
        echo "                              once the installation is complete"
        echo "  -p, --prefix=<path>         installation directory;"
        echo "                              default value is \"../install\""
        echo "  --skip-package-installation don't install third party libraries"
        echo "                              (assume all packages have already been installed"
        echo "                              during a previous execution of this script)"
        echo "  --with-system-<package>     use the default version of <package> installed on"
        echo "                              the system; <package> can be one of:"
        for package in "${PACKAGE_ARRAY[@]}"; do
            echo "                                - $package"
        done
        echo "  --with-system-packages      use the default version installed on the system"
        echo "                              for all packages"
        echo "  --with-<package>=<path>     specify where a package is already installed;"
        echo "                              <package> can be one of:"
        for package in "${PACKAGE_ARRAY[@]}"; do
            echo "                                - $package"
        done
        echo "                              <package> won't be installed; <path> is the"
        echo "                              prefix where the package was installed and"
        echo "                              typically contains the subdirectories bin, lib and"
        echo "                              include (e.g. --with-${PACKAGE_ARRAY[1]}=/usr)"
        echo "  --without-doc               don't build doxygen documentation"
        echo "  --without-tests             don't build tests"
        padded_name=${LOWERCASE_PROJECT_NAME}
        while [ ${#padded_name} -lt 17 ]; do padded_name="${padded_name} "; done
        echo "  --without-${padded_name} only install third party libraries;"
        echo "                              don't install ${PROJECT_NAME}"
        exit 0
    ;;
    --with-system-*)
        PACKAGE_NAME=$(echo $i | sed 's/--with-system-//' | tr '[:upper:]' '[:lower:]')
        if [ "${PACKAGE_NAME}" == "packages" ]; then
            # use the system version for all packages
            for package in "${PACKAGE_ARRAY[@]}"; do
                eval ${package}="SYSTEM" # 'SYSTEM' means the system version of the package will be used
            done
        else
            # use the system version of this package
            eval ${PACKAGE_NAME}="SYSTEM" # 'SYSTEM' means the system version of the package will be used
        fi
    ;;
    --with-*=*)
        PACKAGE_NAME=$(echo $i | sed 's/--with-//' | sed 's/=.*//' | tr '[:upper:]' '[:lower:]')
        eval ${PACKAGE_NAME}="CUSTOM" # 'CUSTOM' means the package is already installed and its path has been specified
        package_install_dir=$(echo $i | sed 's/[-a-zA-Z0-9]*=//')
        # Expand ~ ~+ ~- ~user and remove starting and trailing quotes
        eval package_install_dir=${package_install_dir}
        # If it's not an absolute path prefix it with ${PWD}
        [[ ${package_install_dir} == /* ]] || package_install_dir="${PWD}/${package_install_dir}"
        eval ${PACKAGE_NAME}_install_dir=${package_install_dir}
    ;;
    --without-doc)
        INSTALL_DOC=""
    ;;
    --without-tests)
        BUILD_PROJECT_TESTS="OFF"
    ;;
    --skip-package-installation)
        INSTALL_PACKAGES=""
    ;;
    --without-${LOWERCASE_PROJECT_NAME})
        INSTALL_PROJECT=""
    ;;
    *)
        echo "Error: unrecognized option: $i" >&2
        echo "Try \`$0 --help' for more information." >&2
        exit 1
    ;;
esac
done

################
# Installation #
################

if [ -z "${BUILD_DIR}" ]; then
    BUILD_DIR="${INSTALL_DIR}/build"
fi

# set environment variables to configure compilers
if [ -n "${C_COMPILER}" ]; then
    export CC="${C_COMPILER}"
fi
if [ -n "${CXX_COMPILER}" ]; then
    export CXX="${CXX_COMPILER}"
fi
if [ -n "${FORTRAN_COMPILER}" ]; then
    export FC="${FORTRAN_COMPILER}"
fi

for C_LIBRARY_DIR in "${C_LIBRARY_DIRS[@]}"; do
    export CFLAGS="${CFLAGS} -Wl,-rpath,${C_LIBRARY_DIR} -L${C_LIBRARY_DIR}"
done
for CXX_LIBRARY_DIR in "${CXX_LIBRARY_DIRS[@]}"; do
    export CXXFLAGS="${CXXFLAGS} -Wl,-rpath,${CXX_LIBRARY_DIR} -L${CXX_LIBRARY_DIR}"
done
for FORTRAN_LIBRARY_DIR in "${FORTRAN_LIBRARY_DIRS[@]}"; do
    export FFLAGS="${FFLAGS} -Wl,-rpath,${FORTRAN_LIBRARY_DIR} -L${FORTRAN_LIBRARY_DIR}"
done

echo "-- Build type: ${BUILD_TYPE}"
echo "-- Installing in: $INSTALL_DIR"
echo "-- Using ${NUM_PROC} parallel builds"

# install packages
pushd ${ROOT_DIR}/external > /dev/null
source ./build
export NUM_PROC
popd > /dev/null

CMAKE_BIN=cmake
if [ -n "${cmake_install_dir}" ]; then
    CMAKE_BIN="${cmake_install_dir}/bin/cmake"
fi

if [ -n "${INSTALL_PROJECT}" ]; then
    mkdir -p ${BUILD_DIR}
    cd ${BUILD_DIR}
    echo -n "Configuring ${PROJECT_NAME}... "

    ${CMAKE_BIN} \
        -DTPL_Eigen_DIR:PATH="${eigen_install_dir}" \
        -DTPL_SuiteSparse_DIR:PATH="${suitesparse_install_dir}" \
        -DTPL_BLAS_DIR:PATH="${blas_install_dir}" \
        -DCMAKE_BUILD_TYPE:STRING=${BUILD_TYPE} \
        -DCMAKE_INSTALL_PREFIX:PATH="${INSTALL_DIR}" \
        -D${PROJECT_NAME}_ENABLE_TESTS:BOOL=${BUILD_PROJECT_TESTS} \
        $ROOT_DIR > ${INSTALL_DIR}/${LOWERCASE_PROJECT_NAME}.log 2>&1 || { echo -e "\nError configuring ${PROJECT_NAME}\nCheck log file: ${INSTALL_DIR}/${LOWERCASE_PROJECT_NAME}.log"; exit 1; }
        echo "done."

    if [ -z "${CONFIGURE_ONLY}" ]; then
        echo -n "Building ${PROJECT_NAME}... "
        make -j${NUM_PROC} >> ${INSTALL_DIR}/${LOWERCASE_PROJECT_NAME}.log 2>&1 || { echo -e "\nError building ${PROJECT_NAME}\nCheck log file: ${INSTALL_DIR}/${LOWERCASE_PROJECT_NAME}.log"; exit 1; }
        echo "done."

        if [ -n "${INSTALL_DOC}" ]; then
            echo -n "Creating doxygen documentation... "
            make doc >> ${INSTALL_DIR}/${LOWERCASE_PROJECT_NAME}.log 2>&1 || { echo -e "\nError generating doxygen documentation\nCheck log file: ${INSTALL_DIR}/${LOWERCASE_PROJECT_NAME}.log"; exit 1; }
            echo "done."
        fi

        echo -n "Installing ${PROJECT_NAME}... "
        make -j${NUM_PROC} install >> ${INSTALL_DIR}/${LOWERCASE_PROJECT_NAME}.log 2>&1
        echo "done."

        if [ -z "${LEAVE_BUILD_DIR}" ]; then
            echo -n "Removing build directory... "
            rm -R ${INSTALL_DIR}/build >> ${INSTALL_DIR}/${LOWERCASE_PROJECT_NAME}.log 2>&1
            echo "done."
        fi
    fi
fi

