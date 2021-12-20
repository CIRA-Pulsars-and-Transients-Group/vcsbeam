set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -fPIC -Wall -Wextra")
set(COMPILE_FLAGS "${COMPILE_FLAGS} -DDEBUG")

add_library(vcsbeam STATIC
    ascii_header.c
    performance.c
    filter.c
    form_beam.cu
    jones.c
    buffer.c
    pfb.cu
    calibration.c
    metadata.c)

if (MPI_FOUND AND PAL_FOUND AND PSRFITS_UTILS_FOUND)
    target_sources(vcsbeam PRIVATE beam_psrfits.c)
endif ()

if (PAL_FOUND)
    target_sources(vcsbeam PRIVATE geometry.c)
endif ()

if (VDIFIO_FOUND)
    target_sources(vcsbeam PRIVATE beam_vdif.c)
endif ()

if (HYPERBEAM_FOUND)
    target_sources(vcsbeam PRIVATE primary_beam.c)
endif ()

set_target_properties(vcsbeam
    PROPERTIES PUBLIC_HEADER ${CMAKE_BINARY_DIR}/vcsbeam.h
    COMPILE_FLAGS "${COMPILE_FLAGS}")

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
    cufft)

install(TARGETS vcsbeam
    LIBRARY DESTINATION lib
    PUBLIC_HEADER DESTINATION include)


