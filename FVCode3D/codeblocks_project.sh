#!/bin/bash
# To generate a codeblock project
# Launch it from this directory. Codeblck directory will be ./codeblock/
mkdir codeblock
cd codeblock
cmake .. -G"CodeBlocks - Unix Makefiles"
