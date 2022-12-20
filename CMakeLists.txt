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
find_package(MWALIB REQUIRED)
find_package(MPI)
find_package(PAL)
find_package(CFITSIO)
find_package(PSRFITS_UTILS)
find_package(VDIFIO)
find_package(SIGPROC)
find_package(HYPERBEAM)
find_package(XGPU)

if (PSRFITS_UTILS_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DHAVE_PSRFITS")
endif ()

if (VDIFIO_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DHAVE_VDIFIO")
endif ()

if (SIGPROC_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DHAVE_FILTERBANK")
endif ()

configure_file(include/vcsbeam.h.in vcsbeam.h)

add_subdirectory(src)
add_subdirectory(app)
add_subdirectory(utils)

if(XGPU_FOUND AND CFITSIO_FOUND)
    add_subdirectory(offline_correlator)
endif ()

message("Runtime files will be installed in ${RUNTIME_DIR}")

install(FILES pq_phase_correction.txt
    pfb_filter/MIRROR.dat
    pfb_filter/LSQ12.dat
    pfb_filter/FINEPFB.dat
    DESTINATION ${RUNTIME_DIR})
