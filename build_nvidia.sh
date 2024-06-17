#!/bin/bash

# Dependencies for vcstools
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

ROOT=/software/projects/${PAWSEY_PROJECT}/${USER}
PREFIX=${ROOT}

# PACKAGE is the name of the software
PACKAGE=vcsbeam
VERSION=main_gpu

cd ${ROOT}/${PACKAGE}/src/${PACKAGE}
# Checkout the desired branch
git checkout main_gpu
git pull

# Set the cmake prefix
CMAKE_INSTALL_PREFIX=${ROOT}/modules/${PACKAGE}/${VERSION}



#cd app
#mv make_mwa_incoh_beam.c make_mwa_incoh_beam.cu  
#mv mwa_tied_array_beam_psf.c mwa_tied_array_beam_psf.cu  
#mv mwa_track_primary_beam_response.c mwa_track_primary_beam_response.cu
#mv fine_pfb_offline.c fine_pfb_offline.cu
#cp CMakeLists.txt.NVIDIA CMakeLists.txt
#cd ..

#cd src
#mv pfb.cpp pfb.cu
#mv form_beam.cpp form_beam.cu
#cp CMakeLists.txt.NVIDIA CMakeLists.txt
#cd ..

#cp CMakeLists.txt.NVIDIA CMakeLists.txt

# Build it
mkdir -p build
cd build
#cmake .. -DHYPERBEAM_ROOT=~/github/mwa_hyperbeam/ -DHYPERBEAM_HDF5=/home/msok/github/mwa_hyperbeam/data/mwa_full_embedded_element_pattern.h5 -DVDIFIO_INCLUDE_DIR=/home/msok/github/pulsars/vdifio/src/ -DHYPERBEAM_LIB=/home/msok/github/mwa_hyperbeam/target/release/libmwa_hyperbeam.so -DLIBPAL_LIB=/usr/local/lib/libpal.so
CC=$(which gcc)
CXX=$(which g++)

cmake -DUSE_CUDA=1 -DCMAKE_CUDA_COMPILER=$NVIDIA_CUDA_HOME/bin/nvcc -DCMAKE_CUDA_ARCHITECTURES="native" \
    -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} \
    -DRUNTIME_DIR=${CMAKE_INSTALL_PREFIX}/share/vcsbeam_runtime \
    -DPAL_ROOT_DIR=${PAL_ROOT} \
    -DHYPERBEAM_HDF5=${MWA_BEAM_FILE} \
    -DCFITSIO_ROOT_DIR=${MAALI_CFITSIO_HOME} \
    -DPSRFITS_UTILS_ROOT_DIR=${PSRFITS_UTILS_ROOT} \
    -DVDIFIO_ROOT_DIR=${VDIFIO_ROOT} \
    -DHYPERBEAM_ROOT=${HYPERBEAM_ROOT} \
    -DMWALIB_ROOT=${MWALIB_ROOT} \
    -DXGPU_ROOT=${XGPU_ROOT} \
    ..

make VERBOSE=1
make install

cd ..
rm -rf build
