#! /bin/bash

if [ "$#" -ne 2 ]; then
    echo "$0 <refinement> <input_prefix>"
    exit
fi

TETGEN_DIR="/home/viskius/tetgen1.5.0/build/"
TETGEN_EXE="tetgen"
REFINEMENT=$1
TET_IN="./data/tetgen/in/"$2
TET_IN_FULL=$TET_IN".poly"
TET_OUT="./data/tetgen/out/"

POLYDUAL_DIR="/opt/openfoam231/platforms/linux64GccDPOpt/bin/"
POLYDUAL_EXE="myPolyDualMesh"

# create .node .face and .ele file from .poly file
$TETGEN_DIR/$TETGEN_EXE -pfa$REFINEMENT $TET_IN_FULL
# move output files in the output directory
mv $TET_IN.1.node $TET_IN.1.face $TET_IN.1.ele $TET_OUT
rm -f $TET_IN.1.edge

# convert the tetgen mesh into the openFoam mesh
./polyGen.exe -f tet2foam.txt

pushd .

cd ./data/foam

# create the dual polyhedral mesh using openFoam
$POLYDUAL_DIR/$POLYDUAL_EXE -concaveMultiCells -internalFeatures 0

popd

# convert the openFoam mesh into the fvg mesh
./polyGen.exe -f foam2fvg.txt
