set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -fPIC -Wall -Wextra")
#set(COMPILE_FLAGS "${COMPILE_FLAGS} -DDEBUG")

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

if (PSRFITS_UTILS_FOUND)
    target_compile_options(vcsbeam PRIVATE $<$<COMPILE_LANGUAGE:CUDA>: -DHAVE_PSRFITS>)
endif ()

if (SIGPROC_FOUND)
    target_compile_options(vcsbeam PRIVATE $<$<COMPILE_LANGUAGE:CUDA>: -DHAVE_FILTERBANK>)
endif ()

if (MPI_FOUND AND PAL_FOUND AND PSRFITS_UTILS_FOUND)
    target_sources(vcsbeam PRIVATE beam_psrfits.c)
endif ()

if (PAL_FOUND)
    target_sources(vcsbeam PRIVATE geometry.c)
endif ()

if (VDIFIO_FOUND)
    target_sources(vcsbeam PRIVATE beam_vdif.c)
    target_compile_options(vcsbeam PRIVATE $<$<COMPILE_LANGUAGE:CUDA>: -DHAVE_VDIF>)
endif ()

if (HYPERBEAM_FOUND)
    target_sources(vcsbeam PRIVATE primary_beam.c)
endif ()

set_target_properties(vcsbeam
    PROPERTIES PUBLIC_HEADER ${CMAKE_BINARY_DIR}/vcsbeam.h
    COMPILE_FLAGS "${COMPILE_FLAGS}")

target_include_directories(vcsbeam PUBLIC
    ${CFITSIO_INCLUDE_DIR}
    ${PAL_INCLUDE_DIR}
    ${CUDA_INCLUDE_DIRS}
    ${HYPERBEAM_INCLUDE_DIR}
    ${MWALIB_INCLUDE_DIR}
    ${MPI_INCLUDE_PATH}
    ${CMAKE_BINARY_DIR})

target_link_libraries(vcsbeam
    ${CFITSIO_LIBRARY}
    ${PAL_LIBRARY}
    ${M_LIBRARY}
    ${HYPERBEAM_LIB}
    ${MWALIB_LIB}
    ${MPI_C_LIBRARIES}
    cufft)

if (PSRFITS_UTILS_FOUND)
    target_include_directories(vcsbeam PUBLIC ${PSRFITS_UTILS_INCLUDE_DIR})
    target_link_libraries(vcsbeam ${PSRFITS_UTILS_LIBRARY})
endif ()

if (VDIFIO_FOUND)
    target_include_directories(vcsbeam PUBLIC ${VDIFIO_INCLUDE_DIR})
    target_link_libraries(vcsbeam ${VDIFIO_LIBRARY})
endif ()

if (SIGPROC_FOUND)
    target_include_directories(vcsbeam PUBLIC ${SIGPROC_INCLUDE_DIR})
    target_link_libraries(vcsbeam ${SIGPROC_LIBRARY})
endif ()

install(TARGETS vcsbeam
    LIBRARY DESTINATION lib
    PUBLIC_HEADER DESTINATION include)


