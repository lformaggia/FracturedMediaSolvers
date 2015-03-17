#!/bin/bash
# The archive name can be passed as an argument to this script, e.g.:
# create-archive.sh /path/to/project.tar.gz
# To avoid having to download all packages each time, you can put
# their archives in installer/src

version="HEAD"

root_dir=${PWD}
scripts_dir=installer

source project-config

lowercase_project_name=msrupscaling

output_file="../${lowercase_project_name}.tar.gz"
if [ $# -eq 1 ]; then
    output_file="$1"
fi

# The download_${package} function will look for the archive in
# ${PACKAGE_SRC_DIR} and it will either copy the archive it finds there
# or download it from the internet to ${src_dir}
tmp_dir=tmpdir
src_dir_relative_path=${tmp_dir}/packages
src_dir="${PWD}/${src_dir_relative_path}"
mkdir -p ${src_dir_relative_path}
archive_dir="${scripts_dir}/src"
packages_src_dir="${PWD}/${archive_dir}"
echo "If you interrupt this script before its end, remove the temporary directory:"
echo "  rm -r ${PWD}/${tmp_dir}"
echo "Dowloading or copying packages"

number_of_packages=${#package_array[@]}
for(( i=0; i < ${number_of_packages}; i++ )); do
    package=${package_array[$i]}
    echo "[$((i + 1))/${number_of_packages}] ${package}"
    source ${scripts_dir}/scripts/pkgbuild_${package}
    pushd ${src_dir} > /dev/null
        download_${package}
    popd > /dev/null
done

echo ""
echo -n "Creating an archive of this git repository..."
git archive --format=tar --prefix=${lowercase_project_name}/ -o ${tmp_dir}/project.tar ${version}
echo "done."
echo -n "Appending packages to the archive..."
tar --append -f ${tmp_dir}/project.tar --transform "s,^${src_dir_relative_path},${lowercase_project_name}/${archive_dir}," $(ls -d ${src_dir_relative_path}/*)
echo "done."
echo -n "Compressing the archive..."
gzip -c ${tmp_dir}/project.tar > ${output_file}
echo "done."
echo -n "Removing temporary directory..."
rm -R ${tmp_dir}
echo "done."
