# Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

#------------------------------------------------------------------------------------#
# Install directories
#------------------------------------------------------------------------------------#

# Directory where the CMake configuration file will be installed
set(BITPIT_INSTALL_CMAKE_DIR "lib/cmake/${PROJECT_NAME}-${BITPIT_MAJOR_VERSION}.${BITPIT_MINOR_VERSION}")

#------------------------------------------------------------------------------------#
# Generate BITPITConfig.cmake
#------------------------------------------------------------------------------------#
set(BITPIT_CONFIG_CODE
"set(_bitpit_installed_prefix \"${CMAKE_INSTALL_PREFIX}/${BITPIT_INSTALL_CMAKE_DIR}\")
set(_bitpit_requested_prefix \"\${CMAKE_CURRENT_LIST_DIR}\")
get_filename_component(_bitpit_installed_prefix_full \"\${_bitpit_installed_prefix}\" REALPATH)
get_filename_component(_bitpit_requested_prefix_full \"\${_bitpit_requested_prefix}\" REALPATH)
if (_bitpit_installed_prefix_full STREQUAL _bitpit_requested_prefix_full)
    set(BITPIT_INSTALL_PREFIX \"${CMAKE_INSTALL_PREFIX}\")
else ()
    set(BITPIT_INSTALL_PREFIX \"\${CMAKE_CURRENT_LIST_DIR}\")")

# Construct the proper number of get_filename_component(... PATH)
# calls to compute the installation prefix.
string(REGEX REPLACE "/" ";" _count "${BITPIT_INSTALL_CMAKE_DIR}")
foreach(p ${_count})
    set(BITPIT_CONFIG_CODE "${BITPIT_CONFIG_CODE}
    get_filename_component(BITPIT_INSTALL_PREFIX \"\${BITPIT_INSTALL_PREFIX}\" PATH)")
endforeach()

set(BITPIT_CONFIG_CODE "${BITPIT_CONFIG_CODE}
endif ()")

set(BITPIT_CONFIG_CMAKE_DIR "\${BITPIT_INSTALL_PREFIX}/${BITPIT_INSTALL_CMAKE_DIR}")
get_target_property(BITPIT_DEBUG_POSTFIX ${BITPIT_LIBRARY} DEBUG_POSTFIX)
get_target_property(BITPIT_RELEASE_POSTFIX ${BITPIT_LIBRARY} RELEASE_POSTFIX)

configure_file("BITPITConfig.cmake.in"  "${CMAKE_CURRENT_BINARY_DIR}/BITPITConfig.cmake" @ONLY)

#------------------------------------------------------------------------------------#
# Generate BITPITConfigVersion.cmake.in
#------------------------------------------------------------------------------------#

configure_file("BITPITConfigVersion.cmake.in"  "${CMAKE_CURRENT_BINARY_DIR}/BITPITConfigVersion.cmake" @ONLY)

#------------------------------------------------------------------------------------#
# CMake targets
#------------------------------------------------------------------------------------#
add_custom_target(clean-cmake COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/BITPITConfig.cmake" DESTINATION "${CMAKE_INSTALL_PREFIX}/${BITPIT_INSTALL_CMAKE_DIR}")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/BITPITConfigVersion.cmake" DESTINATION "${CMAKE_INSTALL_PREFIX}/${BITPIT_INSTALL_CMAKE_DIR}")

install(FILES "FindBITPIT.cmake" DESTINATION ${BITPIT_INSTALL_CMAKE_DIR})
install(FILES "UseBITPIT.cmake" DESTINATION "${CMAKE_INSTALL_PREFIX}/${BITPIT_INSTALL_CMAKE_DIR}")
install(FILES "LibFindMacros.cmake" DESTINATION "${CMAKE_INSTALL_PREFIX}/${BITPIT_INSTALL_CMAKE_DIR}")
