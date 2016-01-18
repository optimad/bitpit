# FindBITPIT.cmake
# -----------
#
# The following variables are set if BITPIT is found. If BITPIT is
# not found, BITPIT_FOUND is set to false:
#
#  BITPIT_FOUND        - System has the BITPIT library
#  BITPIT_USE_FILE     - CMake file to use BITPIT.
#  BITPIT_VERSION      - Version of the BITPIT library found
#  BITPIT_INCLUDE_DIRS - The BITPIT include directories
#  BITPIT_LIBRARIES    - The libraries needed to use BITPIT
#  BITPIT_DEFINITIONS  - Compiler switches required for using patchman
#
# The following cache entries must be set by the user to locate BITPIT:
#
#  BITPIT_DIR - The directory containing BITPITConfig.cmake.
#

# Assume not found.
set(BITPIT_FOUND 0)

# Use the Config mode of the find_package() command to find BITPITConfig.
# If this succeeds (possibly because BITPIT_DIR is already set), the
# command will have already loaded BITPITConfig.cmake and set BITPIT_FOUND.
find_package(BITPIT QUIET NO_MODULE)

# If BITPIT was not found, explain to the user how to specify its location.
if (NOT BITPIT_FOUND)
    set(BITPIT_DIR_MESSAGE "BITPIT not found. Set the BITPIT_DIR cmake cache entry to the directory containing BITPITConfig.cmake")

    if (BITPIT_FIND_REQUIRED)
        message(FATAL_ERROR ${BITPIT_DIR_MESSAGE})
    elseif (NOT BITPIT_FIND_QUIETLY)
        message(STATUS ${BITPIT_DIR_MESSAGE})
    endif ()
endif ()
