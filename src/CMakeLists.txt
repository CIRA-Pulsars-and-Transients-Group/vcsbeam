set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -fPIC -Wall -Wextra")

add_library(vcsbeam STATIC
    beam_psrfits.c
    beam_vdif.c
    filter.c
    form_beam.cu
    jones.c
    pfb.cu
    ascii_header.c
    vdifio.c
    calibration.c
    primary_beam.c
    geometry.c
    metadata.c
    performance.c)

set_target_properties(vcsbeam
    PROPERTIES PUBLIC_HEADER ${CMAKE_BINARY_DIR}/vcsbeam.h
    COMPILE_FLAGS "${COMPILE_FLAGS}")

target_include_directories(vcsbeam PUBLIC
    ${PSRFITS_UTILS_INCLUDE_DIR}
    ${CFITSIO_INCLUDE_DIR}
    ${PAL_INCLUDE_DIR}
    ${CUDA_INCLUDE_DIRS}
    ${HYPERBEAM_INCLUDE_DIR}
    ${MWALIB_INCLUDE_DIR}
    ${MPI_INCLUDE_PATH})

target_link_libraries(vcsbeam
    ${PSRFITS_UTILS_LIBRARY}
    ${CFITSIO_LIBRARY}
    ${PAL_LIBRARY}
    ${M_LIBRARY}
    ${HYPERBEAM_LIB}
    ${MWALIB_LIB}
    ${MPI_C_LIBRARIES}
    cufft)

install(TARGETS vcsbeam
    LIBRARY DESTINATION lib
    PUBLIC_HEADER DESTINATION include)


