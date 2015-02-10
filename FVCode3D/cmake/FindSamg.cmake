
if (SAMG_INCLUDES AND SAMG_LIBRARIES)
  set(SAMG_FIND_QUIETLY TRUE)
endif (SAMG_INCLUDES AND SAMG_LIBRARIES)

find_path(SAMG_INCLUDES
  NAMES
  samg.h
  PATHS
  $ENV{SAMGDIR}
  ${INCLUDE_INSTALL_DIR}
  PATH_SUFFIXES
  samg
)

find_library(SAMG_LIBRARIES
  NAMES
  amg
  PATHS
  $ENV{SAMGDIR}
  ${LIB_INSTALL_DIR}
  PATH_SUFFIXES
  samg
)

if(SAMG_LIBRARIES)

  if (NOT SAMG_LIBDIR)
    get_filename_component(SAMG_LIBDIR ${SAMG_LIBRARIES} PATH)
  endif(NOT SAMG_LIBDIR)

  find_library(IOMP5_LIBRARY
    NAMES
    iomp5
    PATHS
    $ENV{SAMGDIR}
    ${LIB_INSTALL_DIR}
    PATH_SUFFIXES
    samg
  )

  if(IOMP5_LIBRARY)
    set(SAMG_LIBRARIES ${SAMG_LIBRARIES} ${IOMP5_LIBRARY})
  else(IOMP5_LIBRARY)
    set(SAMG_LIBRARIES FALSE)
  endif(IOMP5_LIBRARY)

endif(SAMG_LIBRARIES)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SAMG DEFAULT_MSG
                                  SAMG_INCLUDES SAMG_LIBRARIES)

mark_as_advanced(SAMG_INCLUDES SAMG_LIBRARIES)
