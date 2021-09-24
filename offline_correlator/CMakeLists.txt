# Construct the executable
add_executable(offline_correlator offline_correlator.c fourbit.c corr_utils.c)
include_directories(${CFITSIO_INCLUDE_DIR} ${XGPU_INCLUDE_DIRS})
target_link_libraries(offline_correlator ${CFITSIO_LIBRARY} ${M_LIBRARY} ${XGPU_LIBRARY})

# Installation instructions
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/offline_correlator DESTINATION bin)
