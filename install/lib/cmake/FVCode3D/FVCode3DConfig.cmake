# - Config file for FVCode3D
# It defines the following variables
#  FVCode3D_INCLUDE_DIRS - include directories for FVCode3D
#  FVCode3D_LIBRARY_DIRS - library directories for FVCode3D
#  FVCode3D_LIBRARIES    - libraries to link against
#  FVCode3D_TPL_INCLUDE_DIRS - Third party libraries include directories for FVCode3D
#  FVCode3D_TPL_LIBRARY_DIRS - Third party libraries library directories for FVCode3D
#  FVCode3D_TPL_PACKAGES_PROVIDING_USE_FILES - Names of third party libraries with a USE_FILE (which is automatically sourced by this script)
# Additionally, for each PACKAGE in FVCode3D_TPL_PACKAGES_PROVIDING_USE_FILES a variable named TPL_${PACKAGE}_DIR is provided to find the package (often needed to be able to use the corresponding USE_FILE)

# Some packages don't use any variables but rather provide a USE_FILE that sets
# up CMake to link and include the proper files and directories. These packages
# need to pe found again, because their config file sets some variables used by
# their USE_FILE
get_filename_component(FVCode3D_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
include(${FVCode3D_CMAKE_DIR}/FindThirdPartyLibrary.cmake)

set(FVCode3D_TPL_PACKAGES_PROVINDING_USE_FILES )
foreach(PACKAGE ${FVCode3D_TPL_PACKAGES_PROVINDING_USE_FILES})
    find_third_party_library(${PACKAGE})
    include(${PACKAGE}_USE_FILE)
endforeach()

# Compute paths
set(FVCode3D_INCLUDE_DIRS "${FVCode3D_CMAKE_DIR}/../../../include")
set(FVCode3D_LIBRARY_DIRS "${FVCode3D_CMAKE_DIR}/../../")
set(FVCode3D_LIBRARIES fvcode3d;/home/corinne/eniReservoirGITHUB/eniReservoir/install/opt/suitesparse/lib/libumfpack.a;/home/corinne/eniReservoirGITHUB/eniReservoir/install/opt/suitesparse/lib/libcolamd.a;/home/corinne/eniReservoirGITHUB/eniReservoir/install/opt/suitesparse/lib/libamd.a;/home/corinne/eniReservoirGITHUB/eniReservoir/install/opt/suitesparse/lib/libsuitesparseconfig.a;/home/corinne/eniReservoirGITHUB/eniReservoir/install/opt/suitesparse/lib/libcholmod.a;/home/corinne/eniReservoirGITHUB/eniReservoir/install/opt/suitesparse/lib/libamd.a;/home/corinne/eniReservoirGITHUB/eniReservoir/install/opt/suitesparse/lib/libcolamd.a;/home/corinne/eniReservoirGITHUB/eniReservoir/install/opt/suitesparse/lib/libcamd.a;/home/corinne/eniReservoirGITHUB/eniReservoir/install/opt/suitesparse/lib/libccolamd.a;/home/corinne/eniReservoirGITHUB/eniReservoir/install/opt/blas/lib/libblas.so;rt;tet )

set(FVCode3D_TPL_INCLUDE_DIRS "/home/corinne/eniReservoirGITHUB/eniReservoir/install/opt/suitesparse/include;/home/corinne/eniReservoirGITHUB/eniReservoir/install/opt/suitesparse/include;/home/corinne/eniReservoirGITHUB/eniReservoir/install/opt/eigen/include/eigen3")
set(FVCode3D_TPL_LIBRARY_DIRS "")

