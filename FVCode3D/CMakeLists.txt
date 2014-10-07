# Set the minimum required version of cmake for a project
cmake_minimum_required(VERSION 2.8)
# Project name
set(PROJECT_NAME FVCode3D)
project(${PROJECT_NAME})

set (${PROJECT_NAME}_MAJOR_VERSION 0)
set (${PROJECT_NAME}_MINOR_VERSION 1)
set (${PROJECT_NAME}_PATCH_VERSION 0)
set(${PROJECT_NAME}_VERSION
    ${${PROJECT_NAME}_MAJOR_VERSION}.${${PROJECT_NAME}_MINOR_VERSION}.${${PROJECT_NAME}_PATCH_VERSION})
# This is used when configuring files
# (this way we have a variable which is not dependent on the project name)
set(PROJECT_VERSION ${${PROJECT_NAME}_VERSION})

#SET(CMAKE_VERBOSE_MAKEFILE TRUE CACHE BOOL "Verbose Makefile" FORCE)
SET(CMAKE_COLOR_MAKEFILE TRUE CACHE BOOL "Colored Makefile" FORCE)

# This file contains the definition of special functions
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/IncludeFunctions.cmake)

#ENABLE_TESTING()

# compiler stuff
IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
ENDIF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -frounding-math -DSAMG_UNIX_LINUX")

# release and debug mode -Wno-deprecated
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DEIGEN_NO_DEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -pedantic -D_GLIBCXX_DEBUG -O0 -g -gdwarf-3")
set(CMAKE_CXX_FLAGS_PROFILE "-pg")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-Wall -Wextra -pedantic -O2 -DEIGEN_NO_DEBUG")

if("${CMAKE_BUILD_TYPE}" STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release Profile RelWithDebInfo." FORCE)
endif()

#================================
# Installation directories
#================================
# Offer the user the choice of overriding the installation directories
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/../install CACHE PATH "Install directory")
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH
  "Installation directory for header files")

if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR CMake)
else()
  set(DEF_INSTALL_CMAKE_DIR lib/cmake/${PROJECT_NAME})
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH
  "Installation directory for CMake files")

mark_as_advanced(INSTALL_LIB_DIR INSTALL_BIN_DIR INSTALL_INCLUDE_DIR INSTALL_CMAKE_DIR)

#================================
# Third Party Libraries
#================================
# find_package uses a Find<package>.cmake module which it tries to find
# in CMAKE_MODULE_PATH. Since there are no default modules for Cholmod
# and Umfpack, we use those in ${CMAKE_CURRENT_SOURCE_DIR}/cmake
# (we later restore CMAKE_MODULE_PATH its value).
# CMAKE_PREFIX_PATH is where find_package, find_library and find_file
# look for files. We set this to the user defined TPL_<package>_DIR
set(ENV_CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH})
set(ENV_CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH})
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
# CHOLMOD
set(TPL_SuiteSparse_DIR "" CACHE PATH "SuiteSparse installation base directory")
set(CMAKE_PREFIX_PATH ${TPL_SuiteSparse_DIR})
find_package(Cholmod QUIET)
if(NOT CHOLMOD_FOUND)
    message(FATAL_ERROR "Cholmod not found. Try setting TPL_SuiteSparse_DIR to the installation path of SuiteSparse")
endif()
# UMFPACK
set(CMAKE_PREFIX_PATH ${TPL_SuiteSparse_DIR})
find_package(Umfpack QUIET)
if(NOT UMFPACK_FOUND)
    message(FATAL_ERROR "Umfpack not found. Try setting TPL_SuiteSparse_DIR to the installation path of SuiteSparse")
endif()
# Eigen
set(TPL_Eigen_DIR "" CACHE PATH "Eigen installation base directory")
set(CMAKE_PREFIX_PATH ${TPL_Eigen_DIR})
find_package(Eigen3 QUIET)
if(NOT EIGEN3_FOUND)
    message(FATAL_ERROR "Eigen not found. Try setting TPL_Eigen_DIR to the installation path of Eigen")
endif()
# BLAS
set(CMAKE_MODULE_PATH ${ENV_CMAKE_MODULE_PATH})
set(TPL_BLAS_DIR "" CACHE PATH "BLAS installation base directory")
set(CMAKE_PREFIX_PATH ${TPL_BLAS_DIR})
find_package(BLAS QUIET)
if(NOT BLAS_FOUND)
    message(FATAL_ERROR "BLAS not found. Try setting TPL_BLAS_DIR to the installation path of BLAS")
endif()
# Restore variables to their original values
set(CMAKE_PREFIX_PATH ${ENV_CMAKE_PREFIX_PATH})
set(CMAKE_MODULE_PATH ${ENV_CMAKE_MODULE_PATH})

set(${PROJECT_NAME}_TPL_INCLUDE_DIRS
    ${CHOLMOD_INCLUDES}
    ${UMFPACK_INCLUDES}
    ${EIGEN3_INCLUDE_DIR}
    )
message(STATUS "${PROJECT_NAME}_TPL_INCLUDE_DIRS = ${${PROJECT_NAME}_TPL_INCLUDE_DIRS}")

set(${PROJECT_NAME}_TPL_LIBRARY_DIRS
    ""
    )
message(STATUS "${PROJECT_NAME}_TPL_LIBRARY_DIRS = ${${PROJECT_NAME}_TPL_LIBRARY_DIRS}")

# Make sure the compiler can find include files.
include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${${PROJECT_NAME}_TPL_INCLUDE_DIRS}
    )


# Make sure the linker can find the library.
link_directories (
    ${CMAKE_CURRENT_BINARY_DIR}/lib
    ${CMAKE_CURRENT_SOURCE_DIR}/lib
    ${${PROJECT_NAME}_TPL_LIBRARY_DIRS}
    )

set(LIBRARIES
    ${UMFPACK_LIBRARIES}
    ${CHOLMOD_LIBRARIES}
    ${BLAS_LIBRARIES}
    rt
    tet
    samg
    CACHE INTERNAL "")
message(STATUS "LIBRARIES = ${LIBRARIES}")

#================================
# Subdirectories
#================================
option(${PROJECT_NAME}_ENABLE_TESTS "Enable tests" ON)
if(${PROJECT_NAME}_ENABLE_TESTS)
    enable_testing()
endif()
add_subdirectory(src)

# Add documentation
add_subdirectory(doc)

#================================
# ${PROJECT_NAME}Config.cmake
#================================

# This file configures and prepares for installation
# ${PROJECT_NAME}Config.cmake and ${PROJECT_NAME}ConfigVersion.cmake
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/InstallCMakeConfigFiles.cmake)

