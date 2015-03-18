# This function looks for a package in TPL_${PACKAGE_NAME}_DIR first and then
# in the default system paths. If no library is found an error is thrown,
# otherwise its include dirs, library dirs and libraries are added to
# ${PROJECT_NAME}_TPL_INCLUDE_DIRS, ${PROJECT_NAME}_TPL_LIBRARY_DIRS and
# ${PROJECT_NAME}_TPL_LIBRARIES
MACRO(FIND_THIRD_PARTY_LIBRARY PACKAGE_NAME)

# find_package uses a Find<package>.cmake module which it tries to find
# in CMAKE_MODULE_PATH. Since for some packages (umfpack, cholmod, qhull...)
# there is no default module, we use those in ${CMAKE_CURRENT_SOURCE_DIR}/cmake
# (we later restore CMAKE_MODULE_PATH its value).
# CMAKE_PREFIX_PATH is where find_package, find_library and find_file
# look for files. We set this to the user defined TPL_<package>_DIR
# so that system installations of the library don't prevail
set(ENV_CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH})
set(ENV_CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH})
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(TPL_${PACKAGE_NAME}_DIR "" CACHE PATH "${PACKAGE_NAME} installation base directory")
set(CMAKE_PREFIX_PATH ${TPL_${PACKAGE_NAME}_DIR})

# Package specific variables
if("${PACKAGE_NAME}" STREQUAL "Boost")
    set(BOOST_ROOT ${TPL_Boost_DIR})
elseif("${PACKAGE_NAME}" STREQUAL "HDF5")
    set(HDF5_ROOT ${TPL_HDF5_DIR})

elseif("${PACKAGE_NAME}" STREQUAL "ZLIB")
    if(NOT "${TPL_ZLIB_DIR}" STREQUAL "")
        set(ZLIB_ROOT ${TPL_ZLIB_DIR})
    endif("${TPL_ZLIB_DIR}" STREQUAL "")
elseif("${PACKAGE_NAME}" STREQUAL "VTK")
    set(VTK_DIR ${TPL_VTK_DIR})
elseif("${PACKAGE_NAME}" STREQUAL "Trilinos")
    set(Trilinos_DIR ${TPL_Trilinos_DIR}/lib/cmake/Trilinos/)
endif()

# Actually look for the package passing along any extra arguments
find_package(${PACKAGE_NAME} QUIET ${ARGN})

if(NOT ${PACKAGE_NAME}_FOUND)
    message(FATAL_ERROR "${PACKAGE_NAME} not found. Try setting TPL_${PACKAGE_NAME}_DIR to the installation path of ${PACKAGE_NAME}")
endif()

# Restore CMAKE_MODULE_PATH to its original value
set(CMAKE_MODULE_PATH ${ENV_CMAKE_MODULE_PATH})
# Restore CMAKE_PREFIX_PATH to its original value
set(CMAKE_PREFIX_PATH ${ENV_CMAKE_PREFIX_PATH})

# Different packages use different variables to store include directories
if(DEFINED ${PACKAGE_NAME}_INCLUDES)
    set(${PROJECT_NAME}_TPL_INCLUDE_DIRS ${${PROJECT_NAME}_TPL_INCLUDE_DIRS}
        ${${PACKAGE_NAME}_INCLUDES}
        CACHE INTERNAL "")
endif()
if(DEFINED ${PACKAGE_NAME}_INCLUDE_DIR)
    set(${PROJECT_NAME}_TPL_INCLUDE_DIRS ${${PROJECT_NAME}_TPL_INCLUDE_DIRS}
        ${${PACKAGE_NAME}_INCLUDE_DIR}
        CACHE INTERNAL "")
endif()
if(DEFINED ${PACKAGE_NAME}_INCLUDE_DIRS)
    set(${PROJECT_NAME}_TPL_INCLUDE_DIRS ${${PROJECT_NAME}_TPL_INCLUDE_DIRS}
    ${${PACKAGE_NAME}_INCLUDE_DIRS}
    CACHE INTERNAL "")
endif()
if(DEFINED ${PACKAGE_NAME}_INCLUDE_DIR_CPP)
    set(${PROJECT_NAME}_TPL_INCLUDE_DIRS ${${PROJECT_NAME}_TPL_INCLUDE_DIRS}
        ${${PACKAGE_NAME}_INCLUDE_DIR_CPP}
        CACHE INTERNAL "")
endif()

# Library directories
if(DEFINED ${PACKAGE_NAME}_LIBRARY_DIRS)
    set(${PROJECT_NAME}_TPL_LIBRARY_DIRS ${${PROJECT_NAME}_TPL_LIBRARY_DIRS}
        ${${PACKAGE_NAME}_LIBRARY_DIRS}
        CACHE INTERNAL "")
endif()

# Libraries
if(DEFINED ${PACKAGE_NAME}_LIBRARIES)
    set(${PROJECT_NAME}_TPL_LIBRARIES ${${PROJECT_NAME}_TPL_LIBRARIES}
        ${${PACKAGE_NAME}_LIBRARIES}
        CACHE INTERNAL "")
endif()

# Some libraries have a USE_FILE which needs to be included
if(DEFINED ${PACKAGE_NAME}_USE_FILE)
    set(${PROJECT_NAME}_TPL_PACKAGES_PROVIDING_USE_FILES
        ${${PROJECT_NAME}_TPL_PACKAGES_PROVIDING_USE_FILES}
        ${PACKAGE_NAME})
    include(${${PACKAGE_NAME}_USE_FILE})
endif()

# The package itself could have TPL libraries and includes
if(DEFINED ${PACKAGE_NAME}_TPL_INCLUDE_DIRS)
    set(${PROJECT_NAME}_TPL_INCLUDE_DIRS ${${PROJECT_NAME}_TPL_INCLUDE_DIRS}
        ${${PACKAGE_NAME}_TPL_INCLUDE_DIRS}
        CACHE INTERNAL "")
endif()
if(DEFINED ${PACKAGE_NAME}_TPL_LIBRARY_DIRS)
    set(${PROJECT_NAME}_TPL_LIBRAR_DIRS ${${PROJECT_NAME}_TPL_LIBRARY_DIRS}
        ${${PACKAGE_NAME}_TPL_LIBRARY_DIRS}
        CACHE INTERNAL "")
endif()
if(DEFINED ${PACKAGE_NAME}_TPL_PACKAGES_PROVIDING_USE_FILES)
    set(${PROJECT_NAME}_TPL_PACKAGES_PROVIDING_USE_FILES
        ${${PROJECT_NAME}_TPL_PACKAGES_PROVIDING_USE_FILES}
        ${${PACKAGE_NAME}_TPL_PACKAGES_PROVIDING_USE_FILES}
        CACHE INTERNAL "")
endif()

ENDMACRO(FIND_THIRD_PARTY_LIBRARY)
