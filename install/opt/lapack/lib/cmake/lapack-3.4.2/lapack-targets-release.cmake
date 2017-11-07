#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "lapack" for configuration "Release"
set_property(TARGET lapack APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(lapack PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "/home/corinne/eniReservoirGITHUB/eniReservoir/install/opt/blas/lib/libblas.so"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/liblapack.so"
  IMPORTED_SONAME_RELEASE "liblapack.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS lapack )
list(APPEND _IMPORT_CHECK_FILES_FOR_lapack "${_IMPORT_PREFIX}/lib/liblapack.so" )

# Import target "tmglib" for configuration "Release"
set_property(TARGET tmglib APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(tmglib PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "lapack"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libtmglib.so"
  IMPORTED_SONAME_RELEASE "libtmglib.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS tmglib )
list(APPEND _IMPORT_CHECK_FILES_FOR_tmglib "${_IMPORT_PREFIX}/lib/libtmglib.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
