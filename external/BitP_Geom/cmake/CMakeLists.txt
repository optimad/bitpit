# Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

# Directory where to intall the CMake find module
set(DEFAULT_DIR "${BITPIT_DIR}/cmake/Modules")

set(FIND_MODULE_CUSTOMIZE_INSTALL OFF CACHE BOOL "Enables the customization of the directory where the CMake module file will be installed.")
mark_as_advanced(FIND_MODULE_CUSTOMIZE_INSTALL)

if (FIND_MODULE_CUSTOMIZE_INSTALL)
	set(FIND_MODULE_INSTALL_DIR "${DEFAULT_DIR}" CACHE PATH "The directory where the CMake module file will be installed" FORCE)
	mark_as_advanced(FIND_MODULE_INSTALL_DIR)

	set(FIND_MODULE_INSTALL_PATH ${FIND_MODULE_INSTALL_DIR})
else()
	unset(FIND_MODULE_INSTALL_DIR CACHE)
	set(FIND_MODULE_INSTALL_PATH ${DEFAULT_DIR})
endif()

unset(DEFAULT_DIR)

# Generate the CMake find module
string(TOUPPER ${CMAKE_PROJECT_NAME} FIND_MODULE_FILENAME)
set(FIND_MODULE_FILENAME "Find${FIND_MODULE_FILENAME}.cmake")

string(TOUPPER "${BITP_GEOM_LIBRARY}" FIND_MODULE_PACKAGE)
set(FIND_MODULE_NAME "${BITP_GEOM_LIBRARY}")
set(FIND_MODULE_REFERENCE_HEADER "${BITP_GEOM_LIBRARY}_version.hpp")
set(FIND_MODULE_DEPENDENCIES "BITP_MESH")
set(FIND_MODULE_DEPENDENCIES "BITP_BASE")

configure_file("${FIND_MODULE_FILENAME}.in" "${FIND_MODULE_FILENAME}" @ONLY )

# CMake targets
add_custom_target(clean-cmake COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${FIND_MODULE_FILENAME}" DESTINATION "${FIND_MODULE_INSTALL_PATH}")
install(FILES "LibFindMacros.cmake" DESTINATION "${FIND_MODULE_INSTALL_PATH}")
