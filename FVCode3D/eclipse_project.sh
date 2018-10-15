#!/bin/bash
# To generate a codeblock project
# Launch it from this directory. Codeblck directory will be ./codeblock/
cd ../
mkdir eclipse
cd eclipse
cmake ../FVCode3D -G"Eclipse CDT4 - Unix Makefiles"
