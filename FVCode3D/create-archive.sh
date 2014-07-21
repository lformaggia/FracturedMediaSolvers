#!/bin/bash
# The archive name can be passed as an argument to this script, e.g.:
# create-archive.sh /path/to/edfm.tar.gz
# To avoid having to download all packages each time, you can put
# their archives in external/src

LOWERCASE_PROJECT_NAME=fvcode3d

OUTPUT_FILE="../${LOWERCASE_PRJECT_NAME}.tar.gz"
if [ $# -eq 1 ]; then
    OUTPUT_FILE="$1"
fi
PACKAGE_ARRAY=(cmake blas lapack suitesparse eigen)

# The download_${package} function will look for the archive in
# ${PACKAGE_SRC_DIR} and it will either copy the archive it finds there
# or download it from the internet to ${src_dir}
src_dir_relative_path=tmpdir
src_dir="${PWD}/${src_dir_relative_path}"
mkdir -p ${src_dir_relative_path}
archive_dir="external/src"
PACKAGE_SRC_DIR="${PWD}/${archive_dir}"
for package in ${PACKAGE_ARRAY[@]}; do
    source external/scripts/pkgbuild_${package}
    pushd ${src_dir} > /dev/null
        download_${package}
    popd > /dev/null
done
tar -czf ${OUTPUT_FILE} --transform "s,^${src_dir_relative_path},${archive_dir}," --transform "s,^,${LOWERCASE_PROJECT_NAME}/," $(git ls-files) $(ls -d ${src_dir_relative_path}/*)
rm -R ${src_dir_relative_path}
