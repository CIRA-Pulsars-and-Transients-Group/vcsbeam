
# module load rocm/5.4.3
module load rocm/5.4.6
module load cray-hdf5-parallel/1.12.2.7
module load cfitsio/4.3.0
module load rust/1.78.0
module load pal/0.9.8-yyskiux

export LD_LIBRARY_PATH=/scratch/mwavcs/msok/install/mwa_hyperbeam/target/release:/scratch/mwavcs/msok/install/mwalib/target/release:$LD_LIBRARY_PATH

mkdir build
cd build
cmake .. -DUSE_HIP=ON

make VERBOSE=1
