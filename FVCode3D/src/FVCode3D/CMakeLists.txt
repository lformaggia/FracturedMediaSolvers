file(GLOB_RECURSE SRC_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*/*.cpp)
file(GLOB_RECURSE SRC_HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*/*.hpp)
set(SRC_HEADERS ${SRC_HEADERS} FVCode3D.hpp tetgen.h)
# We want these to be cache variables so that they are visible at parent scope
set(LIB_SOURCES ${SRC_SOURCES} CACHE INTERNAL "")
set(LIB_HEADERS ${SRC_HEADERS} CACHE INTERNAL "")

add_library(fvcode3d STATIC ${LIB_SOURCES} ${LIB_HEADERS})

#Set properties on a target.
#We use this here to set -DTETLIBRARY for when compiling the library
set_target_properties(fvcode3d PROPERTIES "COMPILE_DEFINITIONS" TETLIBRARY)

install(TARGETS fvcode3d DESTINATION ${INSTALL_LIB_DIR})
install_files_with_directory_structure(${INSTALL_INCLUDE_DIR}/${PROJECT_NAME} ${LIB_HEADERS})


