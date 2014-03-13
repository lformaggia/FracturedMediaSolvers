#! /bin/bash

# Considering the modules from the Mattia setting
export mkPrefix="/u/geo2/sw"
source ${mkPrefix}/System/mkConfig.sh

# Load the needed modules
module load gcc+system
module load cmake
module load eigen

mkdir -p FVCode3D-build

pushd FVCode3D-build > /dev/null

    rm -rf CMakeCache.txt CMakeFiles

    cmake ../FVCode3D \
        -DEIGEN_PATH:PATH=${mkEigenInc} \
        -DCMAKE_CXX_COMPILER=${mkCxxCompiler} \
        -DCMAKE_C_COMPILER=${mkCCompiler} \
        -DCMAKE_BUILD_TYPE=RELEASE
#        -DCMAKE_BUILD_TYPE=DEBUG

    make

popd > /dev/null
