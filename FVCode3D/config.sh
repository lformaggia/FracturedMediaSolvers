#! /bin/bash
# Usage: ./config.sh RELATIVE_PATH_SOUCE_FOLDER
# To use an absolute path pre-define the variable SOURCE_FOLDER

# Considering the modules from the Mattia setting
export mkPrefix="/u/geo2/sw"
export IS_MOX="true"
source ${mkPrefix}/System/mkConfig.sh

# Load the needed modules
module load gcc+system
module load cmake
module load eigen

# Use a relative path to define the source folder
if [ -z "$SOURCE_FOLDER" ]; then
    SOURCE_FOLDER=${PWD}/${1}
    echo "Usign a relative path: ${SOURCE_FOLDER}"
else
    echo "Usign an absolute path: ${SOURCE_FOLDER}"
fi

NUM_PROC=`nproc`

mkdir -p FVCode3D-build

pushd FVCode3D-build > /dev/null

    rm -rf CMakeCache.txt CMakeFiles

    cmake ${SOURCE_FOLDER} \
        -DEIGEN_PATH:PATH=${mkEigenInc} \
        -DSUITESPARSE_PATH:PATH=${mkSuitesparseHome} \
        -DCMAKE_CXX_COMPILER=${mkCxxCompiler} \
        -DCMAKE_C_COMPILER=${mkCCompiler} \
        -DCMAKE_BUILD_TYPE=DEBUG
#        -DCMAKE_BUILD_TYPE=RELEASE

    make -j${NUM_PROC}

popd > /dev/null

unset SOURCE_FOLDER
