#!/bin/bash
module use /software/projects/pawsey1045/setonix/2024.05/modules/zen3/gcc/12.2.0/
module load mwalib/1.3.3-qvtlpxn
module load pal/0.9.8-yyskiux 
module load cfitsio/4.3.0 psrfits-utils/2023-10-08-ltewgrw vdifio/master-u6heigs  hyperbeam/0.8.0-ub4zrna 
mkdir build
cd build
cmake -DUSE_HIP=ON -DCMAKE_CXX_COMPILER=hipcc -DHYPERBEAM_HDF5=/scratch/mwavcs/msok/install/mwa_hyperbeam/mwa_full_embedded_element_pattern.h5 -DCMAKE_C_COMPILER=hipcc -DPSRFITS_UTILS_ROOT_DIR=${PAWSEY_PSRFITS_UTILS_HOME} -DPAL_ROOT_DIR=${PAWSEY_PAL_HOME} ..
make VERBOSE=1

