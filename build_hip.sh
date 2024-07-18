#!/bin/bash
module use /software/projects/pawsey1045/setonix/2024.05/modules/zen3/gcc/12.2.0/
module load mwalib/1.3.3-qvtlpxn
module load pal/0.9.8-yyskiux 
module load cfitsio/4.3.0 psrfits-utils/2023-10-08-ltewgrw vdifio/master-u6heigs hyperbeam/0.5.0-glmva5q  
mkdir build
cd build
cmake -DUSE_HIP=ON -DCMAKE_INSTALL_PREFIX=/software/projects/pawsey1045/cdipietrantonio/setonix/2024.05/software/manual/vcsbeam/dev -DCMAKE_CXX_COMPILER=hipcc -DHYPERBEAM_HDF5=/software/projects/pawsey1045/setonix/2024.05/software/linux-sles15-zen3/gcc-12.2.0/hyperbeam-0.5.0-glmva5qudgcsa4zv7ufm4nl2id53nqkt/mwa_full_embedded_element_pattern.h5 -DCMAKE_C_COMPILER=hipcc -DPSRFITS_UTILS_ROOT_DIR=${PAWSEY_PSRFITS_UTILS_HOME} -DPAL_ROOT_DIR=${PAWSEY_PAL_HOME} ..
make VERBOSE=1 -j 12
# make install
