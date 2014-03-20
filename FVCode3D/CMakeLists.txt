# Set the minimum required version of cmake for a project
cmake_minimum_required(VERSION 2.8)
# Project name
project (FVCode)			

SET(CMAKE_VERBOSE_MAKEFILE TRUE CACHE BOOL "Verbose Makefile" FORCE)
SET(CMAKE_COLOR_MAKEFILE TRUE CACHE BOOL "Colored Makefile" FORCE)

#ENABLE_TESTING()

# compiler stuff
IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
ENDIF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_FLAGS "-std=c++11 -frounding-math -lumfpack" CACHE STRING "g++ Compiler Flags for All Builds" FORCE)

# release and debug mode -Wno-deprecated
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DEIGEN_NO_DEBUG" CACHE STRING "g++ Compiler Flags for Release" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -pedantic -O0 -g" CACHE STRING "g++ Compiler Flags for Debug" FORCE)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-Wall -Wextra -pedantic -O2 -DEIGEN_NO_DEBUG" CACHE STRING "g++ Compiler Flags for RelWithDebInfo" FORCE)

# install directory
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/../install CACHE PATH "Install directory")
set(INSTALL_BIN ${CMAKE_INSTALL_PREFIX}/bin CACHE PATH "Install: bin subdirectory")
set(INSTALL_LIBRARY ${CMAKE_INSTALL_PREFIX}/lib CACHE PATH "Install: library subdirectory")
set(INSTALL_INCLUDE ${CMAKE_INSTALL_PREFIX}/include CACHE PATH "Install: include subdirectory")

# Eigen path
find_path(EIGEN_PATH NAMES signature_of_eigen3_matrix_library PATHS /usr/local/include/eigen3 DOC "Eigen Directory")

# Suitesparse path
find_path(SUITESPARSE_PATH NAMES signature_of_suitesparse_matrix_library PATHS /usr/local/include/suitesparse DOC "Suitesparse Directory")


# Make sure the compiler can find include files.
include_directories (
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${CMAKE_CURRENT_SOURCE_DIR}/src
            ${CMAKE_CURRENT_SOURCE_DIR}/include
            ${EIGEN_PATH}
            ${SUITESPARSE_PATH}/include
            )


# Make sure the linker can find the library.
link_directories (
            ${CMAKE_CURRENT_BINARY_DIR}/lib
            ${CMAKE_CURRENT_SOURCE_DIR}/lib
            ${SUITESPARSE_PATH}/lib
         )


# Create static library
file(GLOB_RECURSE LIB_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.cpp)
file(GLOB_RECURSE LIB_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.hpp)
add_library(fvcode STATIC ${LIB_SOURCES})

#Set properties on a target. 
#We use this here to set -DTETLIBRARY for when compiling the library
set_target_properties(fvcode PROPERTIES "COMPILE_DEFINITIONS" TETLIBRARY)


# Create executable and link libraries
add_executable(fvcode.exe main.cpp )
target_link_libraries(fvcode.exe fvcode tet)


# Add subdirectories with examples
# add_subdirectory (${CMAKE_CURRENT_SOURCE_DIR}/tests)


# Add documentation
add_subdirectory (${CMAKE_CURRENT_SOURCE_DIR}/doc)

# Add tools
add_subdirectory (${CMAKE_CURRENT_SOURCE_DIR}/tools)

# Install
INSTALL(TARGETS fvcode.exe fvcode
	RUNTIME DESTINATION ${INSTALL_BIN}
	ARCHIVE DESTINATION ${INSTALL_LIBRARY}
	LIBRARY DESTINATION ${INSTALL_LIBRARY}
	)

INSTALL(FILES ${LIB_HEADERS} DESTINATION ${INSTALL_INCLUDE})


# Add data file
add_custom_command(
	TARGET fvcode.exe POST_BUILD        		# Adds a post-build event to fvcode.exe
	
	COMMAND ${CMAKE_COMMAND} -E copy	  	# which executes "cmake - E copy..."
	${PROJECT_SOURCE_DIR}/data.txt     		# <--this is in-file
	${PROJECT_BINARY_DIR}/data.txt			# <--this is out-file path
	
	COMMAND ${CMAKE_COMMAND} -E copy	  	# which executes "cmake - E copy_if_different..."
	${PROJECT_SOURCE_DIR}/README			# <--this is in-file
	${PROJECT_BINARY_DIR}/README			# <--this is out-file path

        COMMAND ln -sf
        "${PROJECT_SOURCE_DIR}/data"
        "${PROJECT_BINARY_DIR}"
	
#	COMMAND ${CMAKE_COMMAND} -E make_directory data

	COMMAND ${CMAKE_COMMAND} -E make_directory results

	)
