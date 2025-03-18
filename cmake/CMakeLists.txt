#---------------------------------------------------------------------------
#
#  bitpit
#
#  Copyright (C) 2015-2023 OPTIMAD engineering Srl
#
#  -------------------------------------------------------------------------
#  License
#  This file is part of bitpit.
#
#  bitpit is free software: you can redistribute it and/or modify it
#  under the terms of the GNU Lesser General Public License v3 (LGPL)
#  as published by the Free Software Foundation.
#
#  bitpit is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
#
#---------------------------------------------------------------------------*/

include(CMakePackageConfigHelpers)

string(TOUPPER ${BITPIT_LIBRARY} UPPER_BITPIT_LIBRARY_NAME)

#------------------------------------------------------------------------------------#
# CMake configuration file
#------------------------------------------------------------------------------------#

set (BITPIT_CMAKE_TARGETS_FILE "${UPPER_BITPIT_LIBRARY_NAME}Targets.cmake")

install(TARGETS ${BITPIT_LIBRARY} EXPORT ${BITPIT_LIBRARY}Targets)

install(EXPORT ${BITPIT_LIBRARY}Targets
    FILE ${BITPIT_CMAKE_TARGETS_FILE}
    NAMESPACE ${BITPIT_LIBRARY}::
    DESTINATION ${BITPIT_INSTALL_CMAKEDIR}
)

set(BITPIT_ENABLED_MODULE_LIST "")
foreach(MODULE_NAME IN LISTS BITPIT_MODULE_LIST)
    isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)
    if(MODULE_ENABLED)
        list(APPEND BITPIT_ENABLED_MODULE_LIST ${MODULE_NAME})
    endif()
endforeach()

set (BITPIT_CMAKE_PACKAGE_FILE "${UPPER_BITPIT_LIBRARY_NAME}Config.cmake")
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/${BITPIT_CMAKE_PACKAGE_FILE}.in
    "${CMAKE_CURRENT_BINARY_DIR}/${BITPIT_CMAKE_PACKAGE_FILE}"
    INSTALL_DESTINATION "${BITPIT_INSTALL_CMAKEDIR}"
    NO_SET_AND_CHECK_MACRO
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

unset(BITPIT_ENABLED_MODULE_LIST)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/${BITPIT_CMAKE_PACKAGE_FILE}
    DESTINATION "${BITPIT_INSTALL_CMAKEDIR}"
)

#------------------------------------------------------------------------------------#
# CMake version file
#------------------------------------------------------------------------------------#

set (BITPIT_CMAKE_CONFIG_VERSION_FILE "${UPPER_BITPIT_LIBRARY_NAME}ConfigVersion.cmake")
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/${BITPIT_CMAKE_CONFIG_VERSION_FILE}
    VERSION ${BITPIT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/${BITPIT_CMAKE_CONFIG_VERSION_FILE}
    DESTINATION "${BITPIT_INSTALL_CMAKEDIR}"
)

#------------------------------------------------------------------------------------#
# CMake find package file
#------------------------------------------------------------------------------------#

install(FILES "FindBITPIT.cmake" DESTINATION "${BITPIT_INSTALL_CMAKEDIR}")

#------------------------------------------------------------------------------------#
# CMake use package file
#------------------------------------------------------------------------------------#

install(FILES "UseBITPIT.cmake" DESTINATION "${BITPIT_INSTALL_CMAKEDIR}")

#------------------------------------------------------------------------------------#
# CMake targets
#------------------------------------------------------------------------------------#
add_custom_target(clean-cmake COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
