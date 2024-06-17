cmake_minimum_required (VERSION 3.15)
#cmake_policy(SET CMP0012 NEW)
#cmake_policy(SET CMP0048 NEW)

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

# Define CMake options to choose the GPU flavour
option(USE_CUDA "Compile the code with NVIDIA GPU support." OFF)
option(USE_HIP "Compile the code with AMD GPU support." OFF)


# Enable CMake support for target language
if(USE_CUDA)
    message("Enabled CUDA (NVIDIA GPU) language support")
    enable_language(CUDA C CXX)
    set(CMAKE_CUDA_ARCHITECTURES "native") 
    add_definitions("-D__NVCC__")
    set_source_files_properties( ${sources} ${headers} ${apps} ${tests} PROPERTIES LANGUAGE CUDA)  # TODO: I don't know if this is correct, what is meant to be in $sources, etc.?
    # this may be needed for older CUDA versions
    add_compile_options("--expt-relaxed-constexpr")
endif()

if(USE_HIP)
    message("Enabled HIP (AMD GPU) language support")
    enable_language(C CXX)
    set(CMAKE_CXX_COMPILER hipcc)
    add_definitions("-D__HIP_PLATFORM_AMD__ -D__HIPCC__")
endif()


project(vcsbeam)

#message("CUDA compute capability is ${SM}")
#if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
#    message("CUDA device compute capability not set -> using default: CMAKE_CUDA_ARCHITECTURES=native")
#    set(CMAKE_CUDA_ARCHITECTURES "native")
#endif()

#set(GPU_LIBS "-lcufft")

#file(GLOB gpu_sources "app/*.cu")
#set_source_files_properties( ${gpu_sources} PROPERTIES LANGUAGE CUDA)

# Add path to hints for package finding
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
include_directories("src/")

# Find packages needed
#find_package(CUDA REQUIRED)
find_package(MWALIB REQUIRED)
find_package(MPI)
find_package(PAL)
find_package(CFITSIO)
find_package(PSRFITS_UTILS)
find_package(HYPERBEAM)
find_package(VDIFIO)
find_package(XGPU)

# Add files/directories to build with
configure_file(include/vcsbeam.h.in vcsbeam.h)
add_subdirectory(src)
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

