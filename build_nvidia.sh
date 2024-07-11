#!/bin/bash
module purge

module load pawseytools/1.29
module load cascadelake/1.0
module load gcc/8.3.0
module load cuda/10.2
module load cfitsio/3450
module load cmake/3.15.0
module use /pawsey/mwa/software/python3/modulefiles
#module load psrfits_utils/master
module load psrfits_utils/284fd0c  # Scott Ransom's version
module load vdifio/master
module load hyperbeam/v0.5.0
module load pal/0.9.8
module load openmpi-ucx-gpu/4.1.6
module load mwalib/v0.16.4
module load xGPU/millisecond


mkdir build_nvidia
cd build_nvidia
cmake -DUSE_CUDA=ON -DHYPERBEAM_HDF5=/pawsey/mwa/mwa_full_embedded_element_pattern.h5 -DPAL_ROOT_DIR=$PAL_ROOT .. 
make VERBOSE=1

