add_executable(fine_pfb_offline fine_pfb_offline.c)
add_executable(make_mwa_tied_array_beam make_mwa_tied_array_beam.c)
add_executable(make_mwa_incoh_beam make_mwa_incoh_beam.c)
add_executable(mwa_track_primary_beam_response mwa_track_primary_beam_response.c)

target_link_libraries(fine_pfb_offline vcsbeam)
target_link_libraries(make_mwa_tied_array_beam vcsbeam)
target_link_libraries(make_mwa_incoh_beam vcsbeam)
target_link_libraries(mwa_track_primary_beam_response vcsbeam)

target_include_directories(fine_pfb_offline PUBLIC ${CMAKE_BINARY_DIR})
target_include_directories(make_mwa_tied_array_beam PUBLIC ${CMAKE_BINARY_DIR})
target_include_directories(make_mwa_incoh_beam PUBLIC ${CMAKE_BINARY_DIR})
target_include_directories(mwa_track_primary_beam_response PUBLIC ${CMAKE_BINARY_DIR})

install(PROGRAMS
    ${CMAKE_CURRENT_BINARY_DIR}/fine_pfb_offline
    ${CMAKE_CURRENT_BINARY_DIR}/make_mwa_tied_array_beam
    ${CMAKE_CURRENT_BINARY_DIR}/make_mwa_incoh_beam
    ${CMAKE_CURRENT_BINARY_DIR}/mwa_track_primary_beam_response
    DESTINATION bin)
