#!/bin/bash
# This script install the development version, meaning that it will not delete
# the build dir and it wil not build the documentation
# Use
#   ./install-devel.sh --help
# for a full list of options

source ./install.sh --leave-build-dir --without-doc $@
