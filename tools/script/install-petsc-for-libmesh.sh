# base directories
PWD=$pwd
local="/home/prashant/Softwares/local_libmesh"
package="/home/prashant/Softwares/Packages"

# For petsc, we simply untar the lib to directory where we want to install library and build files into that directory

lib_name="petsc"
lib_version="3.12.1"

# target info
target_base_dir="$local/$lib_name"

# target build info
debug_build="no"
target_name="opt"
target_dir="$target_base_dir/$lib_version-$target_name"
build_dir="$package/$lib_name/$lib_version-$target_name"

# source info
source_tar_file="$package/$lib_name-$lib_version.tar.gz"
compress_method="gz"

# 
# step 1: Check if various directories exists
#

#
if [[ ! -d $target_base_dir ]]; then 
	mkdir $target_base_dir
fi

#
if [[ ! -d $target_dir ]]; then 
	mkdir $target_dir
fi

#
if [[ ! -d $build_dir ]]; then 
	check_dir="$package/$lib_name"
	if [[ ! -d $check_dir ]]; then
		mkdir $check_dir
	fi
	mkdir $build_dir
else
	cd "$package/$lib_name"
	rm -rf "$lib_version-$target_name"
	mkdir $build_dir
fi

#
# extract tar file
#
tar -zxf "$source_tar_file" -C "$build_dir" --strip-components=1

#
# dependencies
#
metis_dbg_dir="$local/metis/metis-5.1.0/opt"
metis_dir="$local/metis/metis-5.1.0/opt"

# configure and install
cd $build_dir
petsc_arch="arch-linux2-c-debug"
./configure --prefix="$target_dir" \
						--COPTFLAGS='-O3' \
						--CXXOPTFLAGS='-O3' \
						--FOPTFLAGS='-O3' \
						--download-fblaslapack=1 \
						--with-mumps=1 --download-mumps=1 \
						--with-metis=1 --download-metis=1 \
						--with-parmetis=1 --download-parmetis=1 \
						--with-scalapack=1 --download-scalapack=1 \
						--with-superlu=1 --download-superlu=1

make PETSC_DIR="$build_dir" PETSC_ARCH="$petsc_arch" all -j 14

make PETSC_DIR="$build_dir" PETSC_ARCH="$petsc_arch" install

make PETSC_DIR="$target_dir" PETSC_ARCH="" test

make PETSC_DIR="$target_dir" PETSC_ARCH="" streams


# echo info
info_file="$target_dir/install_script.info"
echo "<--- Info --->" | tee -a $info_file
echo "Time = $(date)" | tee -a $info_file

echo ""  | tee -a $info_file
echo "# base directories"  | tee -a $info_file
echo "local = $local" | tee -a $info_file
echo "package = $package" | tee -a $info_file

echo ""  | tee -a $info_file
echo "# library name and version"  | tee -a $info_file
echo "lib_name = $lib_name" | tee -a $info_file
echo "lib_version = $lib_version" | tee -a $info_file
echo "debug_build = $debug_build" | tee -a $info_file
echo "target_name = $target_name" | tee -a $info_file

echo ""  | tee -a $info_file
echo "# target directory, build directory, etc"  | tee -a $info_file
echo "target_base_dir = $target_base_dir" | tee -a $info_file
echo "target_dir = $target_dir" | tee -a $info_file
echo "build_dir = $build_dir" | tee -a $info_file

echo ""  | tee -a $info_file
echo "# source file"  | tee -a $info_file
echo "source_tar_file = $source_tar_file" | tee -a $info_file
echo "compress_method = $compress_method" | tee -a $info_file

echo ""  | tee -a $info_file
echo "# dependencies"  | tee -a $info_file
echo "metis_dbg_dir = $metis_dbg_dir" | tee -a $info_file
echo "metis_dir = $metis_dir" | tee -a $info_file

echo ""  | tee -a $info_file
echo "# PETSC_DIR and PETSC_ARCH"  | tee -a $info_file
echo "PETSC_DIR = $build_dir"  | tee -a $info_file
echo "PETSC_ARCH = $petsc_arch"  | tee -a $info_file

# copy this script to where library is installed
cd $package
cp "install-petsc-for-libmesh.sh" "$target_dir"

