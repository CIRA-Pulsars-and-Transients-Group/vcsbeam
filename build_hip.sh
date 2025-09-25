#!/bin/bash
# TO BUILD ON SETONIX -

# First, you need to source the bash library.

module load bash-utils
source "${BASH_UTILS_DIR}/build_utils.sh"

# Set the program name and versions, used to create the installation paths.
PROGRAM_NAME=vcsbeam
PROGRAM_VERSION="atomics"
# the following function sets up the installation path according to the
# cluster the script is running on and the first argument given. The argument
# can be:
# - "group": install the software in the group wide directory
# - "user": install the software only for the current user
# - "test": install the software in the current working directory
process_build_script_input user

# load all the modules required for the program to compile and run.
# the following command also adds those module names in the modulefile
# that this script will generate.

module_load hyperbeam/0.10.2-cpu mwalib/1.8.7 cfitsio/4.4.0 rocm/6.3.0 pal/0.9.8-gimtz45 psrfits-utils/2023-10-08-bawao2j vdifio/master-xv2az73

module load cmake/3.30.5

mkdir build
cd build
cmake -DUSE_HIP=ON -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR \
    -DCMAKE_CXX_COMPILER=hipcc \
    -DCMAKE_CXX_FLAGS="--offload-arch=gfx90a -O3 -munsafe-fp-atomics" \
    -DHYPERBEAM_HDF5=/scratch/references/mwa/beam-models/mwa_full_embedded_element_pattern.h5 \
    -DCMAKE_C_COMPILER=hipcc \
    -DCMAKE_C_FLAGS="-I${ROCM_PATH}/include" \
    -DCMAKE_BUILD_TYPE=Release \
    -DPSRFITS_UTILS_ROOT_DIR=${PAWSEY_PSRFITS_UTILS_HOME} -DPAL_ROOT_DIR=${PAWSEY_PAL_HOME} ..

make VERBOSE=1 -j 12
make install
create_modulefile
