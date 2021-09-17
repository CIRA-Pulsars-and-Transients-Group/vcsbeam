# - Try to find VDIFIO code.
# Variables used by this module:
#  VDIFIO_ROOT_DIR     - VDIFIO root directory
# Variables defined by this module:
#  VDIFIO_FOUND        - system has VDIFIO
#  VDIFIO_INCLUDE_DIR  - the VDIFIO include directory (cached)
#  VDIFIO_INCLUDE_DIRS - the VDIFIO include directories
#                         (identical to VDIFIO_INCLUDE_DIR)
#  VDIFIO_LIBRARY      - the VDIFIO library (cached)
#  VDIFIO_LIBRARIES    - the VDIFIO libraries

message("Finding VDIFIO")

set(VDIFIO_ROOT_DIR $ENV{VDIFIO})

if(NOT DEFINED VDIFIO_ROOT_DIR)
    message(STATUS "Warning VDIFIO_ROOT_DIR not set: will try and find it ")
else(NOT DEFINED VDIFIO_ROOT_DIR)
    message(STATUS "VDIFIO_ROOT_DIR = ${VDIFIO_ROOT_DIR}")
endif(NOT DEFINED VDIFIO_ROOT_DIR)

if(NOT VDIFIO_FOUND)

  find_path(VDIFIO_INCLUDE_DIR vdifio.h
    HINTS ${VDIFIO_ROOT_DIR} PATH_SUFFIXES /include)
  find_library(VDIFIO_LIBRARY vdifio
    HINTS ${VDIFIO_ROOT_DIR} PATH_SUFFIXES lib )

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(VDIFIO DEFAULT_MSG
    VDIFIO_LIBRARY VDIFIO_INCLUDE_DIR M_LIBRARY)

  set(VDIFIO_INCLUDE_DIRS ${VDIFIO_INCLUDE_DIR})
  set(VDIFIO_LIBRARIES ${VDIFIO_LIBRARY} ${M_LIBRARY})

endif(NOT VDIFIO_FOUND)

if (VDIFIO_FOUND)
    message (STATUS "Found VDIFIO (${VDIFIO_LIBRARIES})")
endif (VDIFIO_FOUND)

