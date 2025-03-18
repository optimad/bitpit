# UseBITPIT.cmake
# -----------
#
# This file sets up include directories, link directories, and
# compiler settings for a project to use BITPIT.  It should not be
# included directly, but rather through the BITPIT_USE_FILE setting
# obtained from BITPITConfig.cmake.
#
# This file is no longer needed and it's here only tom maintain compatibility
# with older versions of BITPIT. It should not be used in new projects.

if(BITPIT_USE_FILE_INCLUDED)
  return()
endif()
set(BITPIT_USE_FILE_INCLUDED 1)

# Warn that the usage of this file is deprecated.
message(WARNING "UseBITPIT.cmake is deprecated and should not be used in new projects.")

# Update CMAKE_MODULE_PATH so includes work.
list(APPEND CMAKE_MODULE_PATH ${BITPIT_CMAKE_DIR})

# Add compiler flags needed to use BITPIT.
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${BITPIT_REQUIRED_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BITPIT_REQUIRED_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${BITPIT_REQUIRED_EXE_LINKER_FLAGS}")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${BITPIT_REQUIRED_SHARED_LINKER_FLAGS}")
set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${BITPIT_REQUIRED_MODULE_LINKER_FLAGS}")

# Add preprocessor definitions needed to use BITPIT.
add_compile_definitions(${BITPIT_DEFINITIONS})

# Add include directories needed to use BITPIT.
include_directories(${BITPIT_INCLUDE_DIRS})
