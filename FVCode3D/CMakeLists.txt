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

SET(CMAKE_VERBOSE_MAKEFILE TRUE CACHE BOOL "Verbose Makefile" FORCE)
SET(CMAKE_COLOR_MAKEFILE TRUE CACHE BOOL "Colored Makefile" FORCE)

# This file contains the definition of special functions
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/IncludeFunctions.cmake)

# compiler stuff
IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
ENDIF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -frounding-math")

# release and debug mode -Wno-deprecated
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DEIGEN_NO_DEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wundef -Wextra -pedantic -D_GLIBCXX_DEBUG -O0 -g -gdwarf-3")
set(CMAKE_CXX_FLAGS_PROFILE "-pg")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-Wall -Wundef -Wextra -pedantic -O2 -DEIGEN_NO_DEBUG")

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
# find_third_party_library looks for a library in TPL_${PACKAGE_NAME}_DIR
# or in the default system paths and appends its include dirs, library dirs
# and library names to the following variables. Some packages don't
# use any variables, but rather provide a USE_FILE that sets up CMake
# to link and include the proper files and directories
set(${PROJECT_NAME}_TPL_INCLUDE_DIRS "" CACHE INTERNAL "")
set(${PROJECT_NAME}_TPL_LIBRARY_DIRS "" CACHE INTERNAL "")
set(${PROJECT_NAME}_TPL_LIBRARIES "" CACHE INTERNAL "")
set(${PROJECT_NAME}_TPL_PACKAGES_PROVIDING_USE_FILES "" CACHE INTERNAL "")

find_third_party_library(UMFPACK)
find_third_party_library(CHOLMOD)
find_third_party_library(BLAS)
find_third_party_library(EIGEN3)
# SAMG
set(TPL_Samg_DIR "" CACHE PATH "Samg installation base directory")
set(CMAKE_PREFIX_PATH ${TPL_Samg_DIR})
find_package(Samg QUIET)
if(NOT SAMG_FOUND)
    message(WARNING "Samg not found. If needed, try setting TPL_Samg_DIR to the installation path of Samg")
    set(SAMG_INCLUDES "")
    set(SAMG_LIBRARIES "")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSAMG_UNIX_LINUX -DSAMG_LCASE_USCORE -DFVCODE3D_HAS_SAMG")
    set(${PROJECT_NAME}_TPL_INCLUDE_DIRS
        ${${PROJECT_NAME}_TPL_INCLUDE_DIRS}
        ${SAMG_INCLUDES}
        CACHE INTERNAL "")
    set(${PROJECT_NAME}_TPL_LIBRARIES
        ${${PROJECT_NAME}_TPL_LIBRARIES}
        ${SAMG_LIBRARIES}
        CACHE INTERNAL "")
endif()

message(STATUS "${PROJECT_NAME}_TPL_INCLUDE_DIRS = ${${PROJECT_NAME}_TPL_INCLUDE_DIRS}")
message(STATUS "${PROJECT_NAME}_TPL_LIBRARY_DIRS = ${${PROJECT_NAME}_TPL_LIBRARY_DIRS}")
message(STATUS "${PROJECT_NAME}_TPL_LIBRARIES = ${${PROJECT_NAME}_TPL_LIBRARIES}")
message(STATUS "${PROJECT_NAME}_TPL_PACKAGES_PROVIDING_USE_FILES = ${${PROJECT_NAME}_TPL_PACKAGES_PROVIDING_USE_FILES}")

include_directories(${${PROJECT_NAME}_TPL_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/include)

link_directories(${${PROJECT_NAME}_TPL_LIBRARY_DIRS}
    ${CMAKE_CURRENT_BINARY_DIR}/lib
    ${CMAKE_CURRENT_SOURCE_DIR}/lib)

set(LIBRARIES
    ${${PROJECT_NAME}_TPL_LIBRARIES}
    rt
    tet
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

