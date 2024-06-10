#!/bin/bash

cd app
mv make_mwa_incoh_beam.c make_mwa_incoh_beam.cu  
mv mwa_tied_array_beam_psf.c mwa_tied_array_beam_psf.cu  
mv mwa_track_primary_beam_response.c mwa_track_primary_beam_response.cu
cp CMakeLists.txt.NVIDIA CMakeLists.txt

cd ../src/
mv pfb.cpp pfb.cu
mv form_beam.cpp form_beam.cu
cp CMakeLists.txt.NVIDIA CMakeLists.txt
cd ..

cp CMakeLists.txt.NVIDIA CMakeLists.txt
mkdir build
cd build
cmake .. -DHYPERBEAM_ROOT=~/github/mwa_hyperbeam/ -DHYPERBEAM_HDF5=/home/msok/github/mwa_hyperbeam/data/mwa_full_embedded_element_pattern.h5 -DVDIFIO_INCLUDE_DIR=/home/msok/github/pulsars/vdifio/src/ -DHYPERBEAM_LIB=/home/msok/github/mwa_hyperbeam/target/release/libmwa_hyperbeam.so -DLIBPAL_LIB=/usr/local/lib/libpal.so

make VERBOSE=1

