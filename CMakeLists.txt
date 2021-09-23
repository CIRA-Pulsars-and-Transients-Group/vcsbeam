cmake_minimum_required (VERSION 3.15)
#cmake_policy(SET CMP0012 NEW)
#cmake_policy(SET CMP0048 NEW)

# Set up version number
execute_process(
  COMMAND bash -c "git describe --tags --long | sed 's/-/./' | sed 's/-g/_/'"
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE VCSBEAM_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

message("VCSBeam ${VCSBEAM_VERSION}")

if(NOT RUNTIME_DIR)
    set(RUNTIME_DIR ${CMAKE_INSTALL_PREFIX}/bin/vcsbeam_runtime)
endif(NOT RUNTIME_DIR)

project (vcsbeam LANGUAGES C CUDA)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(CUDA REQUIRED)
find_package(MPI REQUIRED)
find_package(PAL REQUIRED)
find_package(CFITSIO REQUIRED)
find_package(PSRFITS_UTILS REQUIRED)
find_package(HYPERBEAM REQUIRED)
find_package(MWALIB REQUIRED)
find_package(VDIFIO REQUIRED)
find_package(XGPU)

configure_file(include/vcsbeam.h.in vcsbeam.h)

add_subdirectory(src)
add_subdirectory(app)

message("Runtime files will be installed in ${RUNTIME_DIR}")

install(FILES pq_phase_correction.txt
    pfb_filter/MIRROR.dat
    pfb_filter/LSQ12.dat
    pfb_filter/FINEPFB.dat
    DESTINATION ${RUNTIME_DIR})
