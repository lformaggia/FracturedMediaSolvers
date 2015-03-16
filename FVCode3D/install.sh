#!/bin/bash
# This script installs this project and all the needed dependencies
# For a simple installation:
#   ./install.sh --prefix=<install dir>
# For a full list of option:
#   ./install.sh --help

source project-config

lowercase_project_name=$(echo ${project_name} | tr '[:upper:]' '[:lower:]')
uppercase_project_name=$(echo ${project_name} | tr '[:lower:]' '[:upper:]')

###################
# Default options #
###################

root_dir=$PWD
scripts_dir=${root_dir}/installer

# Basic options
pushd ${scripts_dir} > /dev/null
source config
set_default_compiler_variables
popd > /dev/null

# Project specific options
build_type="Release"
install_project=1
install_doc=1
build_project_tests="ON"
leave_project_build_dir=""

####################
# Argument parsing #
####################
source ${scripts_dir}/argument-parsing
for i in "$@"
    do
        if handle_basic_argument $i; then
            :
        else
            case $i in
                -h|--help)
                    print_basic_help
                    echo ""
                    echo "Build options:"
                    echo "  --build-type=<type>         build ${project_name} in <type> mode,"
                    echo "                              where <type> can be one of:"
                    echo "                               - Release (default)"
                    echo "                               - Debug"
                    echo "                               - Profile"
                    echo "  --configure-only            configure ${project_name} and create makefiles"
                    echo "                              in the build directory, but don't execute"
                    echo "                              the actual build"
                    echo "  --leave-build-dir           ${project_name} build directory won't be deleted"
                    echo "                              once the installation is complete"
                    echo "  --without-doc               don't build doxygen documentation"
                    echo "  --without-tests             don't build tests"
                    padded_name=${lowercase_project_name}
                    while [ ${#padded_name} -lt 17 ]; do padded_name="${padded_name} "; done
                    echo "  --without-${padded_name} only install third party libraries;"
                    echo "                              don't install ${project_name}"
                    echo ""
                    echo "  -h, --help                  print this help"
                    exit 0
                ;;
                --build-type=*)
                    arg_build_type=$(echo $i | sed 's/[-a-zA-Z0-9]*=//')
                    if [ "${arg_build_type,,}" != "release" ] && [ "${arg_build_type,,}" != "debug" ] && [ "${arg_build_type,,}" != "profile" ]; then
                        echo "Unrecognized build type."
                        echo "Valid options are: Release, Debug and Profile"
                        exit 1
                    fi
                    build_type=${arg_build_type}
                ;;
                --configure-only)
                    configure_only=1
                ;;
                --leave-build-dir)
                    leave_project_build_dir=1
                ;;
                --without-doc)
                    install_doc=""
                ;;
                --without-tests)
                    build_project_tests="OFF"
                ;;
                --without-${lowercase_project_name})
                    install_project=""
                ;;
                *)
                    echo "Error: unrecognized option: $i" >&2
                    echo "Try \`$0 --help' for more information." >&2
                    exit 1
                ;;
            esac
        fi
done

# Set the default values for some basic options, based on the argument just parsed
pushd ${scripts_dir} > /dev/null
set_default_basic_options
popd > /dev/null

################
# Installation #
################

# Set environment variables to configure compilers
export_compiler_variables

echo "-- Build type: ${build_type}"
echo "-- Installing in: ${project_install_dir}"
echo "-- Using ${num_proc} parallel builds"

# Install packages
pushd ${scripts_dir} > /dev/null
source ./build
export num_proc
popd > /dev/null

# Install this project
log_file=${project_build_dir}/${lowercase_project_name}.log

cmake_bin=cmake
if [ -n "${cmake_install_dir}" ]; then
    cmake_bin="${cmake_install_dir}/bin/cmake"
fi

if [ -n "${install_project}" ]; then
    mkdir -p ${project_build_dir}
    cd ${project_build_dir}
    echo ""
    echo -n "Configuring ${project_name}... "

    ${cmake_bin} \
        -DTPL_BLAS_DIR:PATH="${blas_install_dir}" \
        -DTPL_EIGEN3_DIR:PATH="${eigen_install_dir}" \
        -DTPL_CHOLMOD_DIR:PATH="${suitesparse_install_dir}" \
        -DTPL_UUMFPACK_DIR:PATH="${suitesparse_install_dir}" \
        -DCMAKE_build_type:STRING=${build_type} \
        -DCMAKE_INSTALL_PREFIX:PATH="${project_install_dir}" \
        -D${project_name}_ENABLE_TESTS:BOOL=${build_project_tests} \
        $root_dir > ${log_file} 2>&1 || { echo -e "\nError configuring ${project_name}\nCheck log file: ${log_file}"; exit 1; }
        echo "done."

    if [ -z "${configure_only}" ]; then
        echo -n "Building ${project_name}... "
        make -j${num_proc} >> ${log_file} 2>&1 || { echo -e "\nError building ${project_name}\nCheck log file: ${log_file}"; exit 1; }
        echo "done."

        if [ -n "${install_doc}" ]; then
            echo -n "Creating doxygen documentation... "
            make doc >> ${log_file} 2>&1 || { echo -e "\nError generating doxygen documentation\nCheck log file: ${log_file}"; exit 1; }
            echo "done."
        fi

        echo -n "Installing ${project_name}... "
        make -j${num_proc} install >> ${log_file} 2>&1 || { echo -e "\nError installing ${project_name}\nCheck log file: ${log_file}"; exit 1; }

        echo "done."

        if [ -z "${leave_project_build_dir}" ]; then
            echo -n "Removing build directory... "
            rm -R ${project_install_dir}/build >> ${log_file} 2>&1
            echo "done."
        fi
    fi
fi

