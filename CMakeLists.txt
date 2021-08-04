cmake_minimum_required (VERSION 3.15)
#cmake_policy(SET CMP0012 NEW)
#cmake_policy(SET CMP0048 NEW)

project (vcsbeam LANGUAGES C CUDA)

# Set up version number
execute_process(
  COMMAND bash -c "git describe --tags --long | sed 's/-/./' | sed 's/-/_/'"
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE VCSBEAM_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

message("VCSBeam ${VCSBEAM_VERSION}")

add_definitions( -DVERSION_BEAMFORMER="${VCSBEAM_VERSION}" )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(CUDA)
find_package(PAL)
find_package(CFITSIO)
find_package(PSRFITS_UTILS)
#find_package(FFTW3 COMPONENTS single openmp)
find_package(HYPERBEAM)


#find_package(OpenMP)
#if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#else()
#    message("No OpenMP")
#endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -fPIC -Wall -Wextra")

set(sources "beam_common.c;beam_psrfits.c;beam_vdif.c;filter.c;form_beam.cu;get_delays_small.c;ipfb.cu;make_beam.c;ascii_header.c;vdifio.c")

add_executable(make_beam ${sources})

set_target_properties(make_beam PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DHAVE_CUDA -DHYPERBEAM_HDF5=\\\"${HYPERBEAM_HDF5}\\\"")

include_directories(${PSRFITS_UTILS_INCLUDE_DIR} ${CFITSIO_INCLUDE_DIR} ${PAL_INCLUDE_DIR} ${CUDA_INCLUDE_DIRS} ${HYPERBEAM_INCLUDE_DIR})

target_link_libraries(make_beam ${PSRFITS_UTILS_LIBRARY} ${CFITSIO_LIBRARY} ${PAL_LIBRARY} ${M_LIBRARY} ${HYPERBEAM_LIB})

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/make_beam DESTINATION bin)

