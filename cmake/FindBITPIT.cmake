# FindBITPIT.cmake
# -----------
#
# The following variables are set if bitpit is found. If bitpit is
# not found, BITPIT_FOUND is set to false:
#
#  BITPIT_FOUND        - System has the bitpit library
#  BITPIT_USE_FILE     - CMake file to use bitpit.
#  BITPIT_VERSION      - Version of the bitpit library found
#  BITPIT_INCLUDE_DIRS - The bitpit include directories
#  BITPIT_LIBRARIES    - The libraries needed to use bitpit
#  BITPIT_DEFINITIONS  - Compiler switches required for using bitpit
#  BITPIT_LANGUAGES    - Programming languages used by the bitpit library
#
# The imported target allows external projects to use the bitpit library
# by simply calling target_link_libraries:
#
#     find_package(bitpit)
#     target_link_libraries(external_target PRIVATE bitpit::bitpit)
#
# The following cache entries must be set by the user to locate bitpit:
#
#  BITPIT_DIR - The directory containing BITPITConfig.cmake.
#
# A list of required bitpit modules may be specified when invoking the
# find_package command after the COMPONENTS option (or after the REQUIRED
# option if present). Additional optional components may be listed after
# OPTIONAL_COMPONENTS. For each of the requested modules, a boolean variable
# named BITPIT_<MODULE_NAME>_FOUND will be set telling if the corresponding
# module is enabled in the corresponding bitpit installation. If a required
# module is not found a fatal error is generated and the configure step
# stops executing.

# Assume not found.
set(BITPIT_FOUND 0)

# Warn that the usage of this file is deprecated.
message(WARNING "FindBITPIT.cmake is deprecated and should not be used in new projects.")

# Use the Config mode of the find_package() command to find BITPITConfig.
# If this succeeds (possibly because BITPIT_DIR is already set), the
# command will have already loaded BITPITConfig.cmake and set BITPIT_FOUND.
find_package(BITPIT QUIET NO_MODULE COMPONENTS ${BITPIT_FIND_COMPONENTS} OPTIONAL_COMPONENTS ${BITPIT_FIND_OPTIONAL_COMPONENTS})

# If bitpit was not found, explain to the user how to specify its location.
if (NOT BITPIT_FOUND)
    set(BITPIT_DIR_MESSAGE "BITPIT not found. Set the BITPIT_DIR cmake cache entry to the directory containing BITPITConfig.cmake")

    if (BITPIT_FIND_REQUIRED)
        message(FATAL_ERROR ${BITPIT_DIR_MESSAGE})
    elseif (NOT BITPIT_FIND_QUIETLY)
        message(STATUS ${BITPIT_DIR_MESSAGE})
    endif ()
endif ()

# If a required module is not found a fatal error is generated and the
# configure step stops executing.
foreach(COMPONENT ${BITPIT_FIND_COMPONENTS})
    if(NOT BITPIT_${COMPONENT}_FOUND)
        set(COMPONENT_NOT_FOUND_MESSAGE "${COMPONENT} module is not enabled in current bitpit installation")
        if(BITPIT_FIND_REQUIRED_${COMPONENT})
           message(FATAL_ERROR "${COMPONENT_NOT_FOUND_MESSAGE}")
        elseif (NOT BITPIT_FIND_QUIETLY)
           message(STATUS "${COMPONENT_NOT_FOUND_MESSAGE}")
        endif ()
    endif()
endforeach()

# If an optional module is not found a simple warning is generated
if (NOT BITPIT_FIND_QUIETLY)
    foreach(COMPONENT ${BITPIT_FIND_OPTIONAL_COMPONENTS})
        if(NOT BITPIT_${COMPONENT}_FOUND)
            set(COMPONENT_NOT_FOUND_MESSAGE "${COMPONENT} optional module is not enabled in current bitpit installation")
            message(STATUS "${COMPONENT_NOT_FOUND_MESSAGE}")
        endif ()
    endforeach()
endif()
