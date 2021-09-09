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

add_definitions( -DVERSION_BEAMFORMER="${VCSBEAM_VERSION}" -DRUNTIME_DIR="${RUNTIME_DIR}" )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(CUDA REQUIRED)
find_package(MPI REQUIRED)
find_package(PAL REQUIRED)
find_package(CFITSIO REQUIRED)
find_package(PSRFITS_UTILS REQUIRED)
find_package(HYPERBEAM REQUIRED)
find_package(MWALIB REQUIRED)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -fPIC -Wall -Wextra")

add_executable(make_mwa_tied_array_beam "beam_psrfits.c;beam_vdif.c;filter.c;form_beam.cu;jones.c;pfb.cu;make_mwa_tied_array_beam.c;ascii_header.c;vdifio.c;calibration.c;primary_beam.c;geometry.c;metadata.c;performance.c")
add_executable(make_mwa_incoh_beam "beam_psrfits.c;form_beam.cu;make_mwa_incoh_beam.c;ascii_header.c;metadata.c;geometry.c;performance.c")
add_executable(mwa_track_primary_beam_response "mwa_track_primary_beam_response.c;jones.c;metadata.c;geometry.c;primary_beam.c;form_beam.cu;performance.c")
add_executable(fine_pfb_offline "fine_pfb_offline.c;metadata.c;performance.c;pfb.cu")

set_target_properties(make_mwa_tied_array_beam PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DHYPERBEAM_HDF5=\\\"${HYPERBEAM_HDF5}\\\"")
set_target_properties(make_mwa_incoh_beam PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
set_target_properties(mwa_track_primary_beam_response PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DHYPERBEAM_HDF5=\\\"${HYPERBEAM_HDF5}\\\"")
set_target_properties(fine_pfb_offline PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")

include_directories(${PSRFITS_UTILS_INCLUDE_DIR} ${CFITSIO_INCLUDE_DIR} ${PAL_INCLUDE_DIR} ${CUDA_INCLUDE_DIRS} ${HYPERBEAM_INCLUDE_DIR} ${MWALIB_INCLUDE_DIR} ${MPI_INCLUDE_PATH})

target_link_libraries(make_mwa_tied_array_beam ${PSRFITS_UTILS_LIBRARY} ${CFITSIO_LIBRARY} ${PAL_LIBRARY} ${M_LIBRARY} ${HYPERBEAM_LIB} ${MWALIB_LIB} ${MPI_C_LIBRARIES} cufft)
target_link_libraries(make_mwa_incoh_beam ${PSRFITS_UTILS_LIBRARY} ${CFITSIO_LIBRARY} ${PAL_LIBRARY} ${M_LIBRARY} ${MWALIB_LIB} ${MPI_C_LIBRARIES})
target_link_libraries(mwa_track_primary_beam_response ${PAL_LIBRARY} ${MWALIB_LIB} ${HYPERBEAM_LIB} ${M_LIBRARY})
target_link_libraries(fine_pfb_offline ${MWALIB_LIB} ${M_LIBRARY} cufft)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/make_mwa_tied_array_beam ${CMAKE_CURRENT_BINARY_DIR}/make_mwa_incoh_beam ${CMAKE_CURRENT_BINARY_DIR}/mwa_track_primary_beam_response DESTINATION bin)

install(FILES xy_phase_correction.txt
    pfb_filter/MIRROR.dat
    pfb_filter/LSQ12.dat
    pfb_filter/FINEPFB.dat
    DESTINATION ${RUNTIME_DIR})

message("Runtime files will be installed in ${RUNTIME_DIR}")
