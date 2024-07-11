cmake_minimum_required (VERSION 3.15)
project(vcsbeam)
#cmake_policy(SET CMP0012 NEW)
#cmake_policy(SET CMP0048 NEW)

# Define CMake options to choose the GPU flavour
option(USE_CUDA "Compile the code with NVIDIA GPU support." OFF)
option(USE_HIP "Compile the code with AMD GPU support." OFF)

# Find packages needed
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
#find_package(CUDA REQUIRED)
find_package(MWALIB REQUIRED)
find_package(MPI REQUIRED)
find_package(PAL REQUIRED)
find_package(CFITSIO REQUIRED)
find_package(PSRFITS_UTILS REQUIRED)
find_package(HYPERBEAM REQUIRED)
find_package(VDIFIO REQUIRED)
find_package(XGPU)

if(USE_CUDA)
    enable_language(CUDA C CXX)
    set(CMAKE_CUDA_ARCHITECTURES "native") 
    add_definitions("-D__NVCC__")
    set(GPU_FFTLIB cufft)
elseif(USE_HIP)
    #  only available from CMAKE 3.21. For versions less than 3.21, comment out the
    # enable_language macro and simply specify the -DCMAKE_CXX_COMPILER=hipcc 
    # command line option.
    # enable_language(HIP CXX)
    set(GPU_FFTLIB hipfft)
    add_definitions("-D__HIP_PLATFORM_AMD__ -D__HIPCC__")
endif()

# Set up version number from Git
execute_process(
  COMMAND bash -c "git describe --tags --long | sed 's/-/./' | sed 's/-g/_/'"
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE VCSBEAM_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
message("VCSBeam ${VCSBEAM_VERSION}")


# Ensure that runtime files have an installation path
if(NOT RUNTIME_DIR)
    set(RUNTIME_DIR ${CMAKE_INSTALL_PREFIX}/bin/vcsbeam_runtime)
endif()
message("Runtime files will be installed in ${RUNTIME_DIR}")

file(GLOB vcsbeam_c_sources 
    "src/ascii_header.c"
    "src/performance.c"
    "src/filter.c"
    "src/jones.c"
    "src/buffer.c"
    "src/calibration.c"
    "src/metadata.c"
)

file(GLOB vcsbeam_gpu_sources 
    "src/form_beam.cpp"
    "src/pfb.cpp"
)

# Enable CMake support for target language
if(USE_CUDA)
    message("Enabled CUDA (NVIDIA GPU) language support")
    set_source_files_properties(${vcsbeam_gpu_sources} PROPERTIES LANGUAGE CUDA)
endif()

add_library(vcsbeam STATIC
    ${vcsbeam_c_sources}
    ${vcsbeam_gpu_sources}
)

if (MPI_FOUND AND PAL_FOUND AND PSRFITS_UTILS_FOUND)
    target_sources(vcsbeam PRIVATE "src/beam_psrfits.c")
endif ()

if (PAL_FOUND)
    target_sources(vcsbeam PRIVATE "src/geometry.c")
endif ()

if (VDIFIO_FOUND)
    target_sources(vcsbeam PRIVATE "src/beam_vdif.c")
endif ()

if (HYPERBEAM_FOUND)
    target_sources(vcsbeam PRIVATE "src/primary_beam.c")
endif ()


target_include_directories(vcsbeam PUBLIC
    ${PSRFITS_UTILS_INCLUDE_DIR}
    ${VDIFIO_INCLUDE_DIR}
    ${CFITSIO_INCLUDE_DIR}
    ${PAL_INCLUDE_DIR}
    ${CUDA_INCLUDE_DIRS}
    ${HYPERBEAM_INCLUDE_DIR}
    ${MWALIB_INCLUDE_DIR}
    ${MPI_INCLUDE_PATH}
    ${CMAKE_BINARY_DIR})

target_link_libraries(vcsbeam
    ${PSRFITS_UTILS_LIBRARY}
    ${VDIFIO_LIBRARY}
    ${CFITSIO_LIBRARY}
    ${PAL_LIBRARY}
    ${M_LIBRARY}
    ${HYPERBEAM_LIB}
    ${MWALIB_LIB}
    ${MPI_C_LIBRARIES}
    ${GPU_FFTLIB})

install(TARGETS vcsbeam
    LIBRARY DESTINATION lib
    PUBLIC_HEADER DESTINATION include)


# Add files/directories to build with
configure_file(include/vcsbeam.h.in ${CMAKE_BINARY_DIR}/vcsbeam.h)

set_target_properties(vcsbeam
    PROPERTIES PUBLIC_HEADER "vcsbeam.h"
)

# Add path to hints for package finding

include_directories("src/")


add_subdirectory(app)
add_subdirectory(utils)

if(XGPU_FOUND AND CFITSIO_FOUND)
    add_subdirectory(offline_correlator)
endif ()


# Install any necessary data files in the required location
install(FILES pq_phase_correction.txt
    pfb_filter/MIRROR.dat
    pfb_filter/LSQ12.dat
    pfb_filter/FINEPFB.dat
    DESTINATION ${RUNTIME_DIR})

