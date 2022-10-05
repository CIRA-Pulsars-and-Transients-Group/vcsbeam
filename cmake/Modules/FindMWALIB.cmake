# - Try to find MWALIB code.
# Variables used by this module:
#  MWALIB_ROOT     - MWALIB root directory
# Variables defined by this module:
#  MWALIB_FOUND        - system has MWALIB
#  MWALIB_INCLUDE_DIR  - the MWALIB include directory (cached)
#  MWALIB_LIB          - the MWALIB library (cached)

message("Finding MWALIB")

set(MWALIB_ROOT $ENV{MWALIB})

if(NOT DEFINED MWALIB_ROOT)
    message(STATUS "Warning MWALIB_ROOT not set: will try and find it ")
else(NOT DEFINED MWALIB_ROOT)
    message(STATUS "MWALIB_ROOT = ${MWALIB_ROOT}")
endif(NOT DEFINED MWALIB_ROOT)

if(NOT MWALIB_FOUND)

  find_path(MWALIB_INCLUDE_DIR mwalib.h
    HINTS ${MWALIB_ROOT} PATH_SUFFIXES /include)
  find_library(MWALIB_LIB mwalib
    HINTS ${MWALIB_ROOT} PATH_SUFFIXES /target/release )

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(MWALIB DEFAULT_MSG
    MWALIB_LIB MWALIB_INCLUDE_DIR)

endif(NOT MWALIB_FOUND)

if (MWALIB_FOUND)
    message (STATUS "Found MWALIB (${MWALIB_LIB})")
endif (MWALIB_FOUND)

