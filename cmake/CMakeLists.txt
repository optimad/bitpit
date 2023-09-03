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

#------------------------------------------------------------------------------------#
# Generate BITPITConfig.cmake
#------------------------------------------------------------------------------------#
set(BITPIT_INSTALL_PREFIX_CODE
"set(_bitpit_installed_prefix \"${BITPIT_INSTALL_CMAKEDIR}\")
set(_bitpit_requested_prefix \"\${CMAKE_CURRENT_LIST_DIR}\")
get_filename_component(_bitpit_installed_prefix_full \"\${_bitpit_installed_prefix}\" REALPATH)
get_filename_component(_bitpit_requested_prefix_full \"\${_bitpit_requested_prefix}\" REALPATH)
if (_bitpit_installed_prefix_full STREQUAL _bitpit_requested_prefix_full)
    set(BITPIT_INSTALL_PREFIX \"${CMAKE_INSTALL_PREFIX}\")
else ()
    set(BITPIT_INSTALL_PREFIX \"\${CMAKE_CURRENT_LIST_DIR}\")")

# Construct the proper number of get_filename_component(... PATH)
# calls to compute the installation prefix.
string(REGEX REPLACE "/" ";" _count "${BITPIT_INSTALL_CMAKEDIR}")
foreach(p ${_count})
    set(BITPIT_INSTALL_PREFIX_CODE "${BITPIT_INSTALL_PREFIX_CODE}
    get_filename_component(BITPIT_INSTALL_PREFIX \"\${BITPIT_INSTALL_PREFIX}\" PATH)")
endforeach()

set(BITPIT_INSTALL_PREFIX_CODE "${BITPIT_INSTALL_PREFIX_CODE}
endif ()")

get_target_property(BITPIT_DEBUG_POSTFIX ${BITPIT_LIBRARY} DEBUG_POSTFIX)
get_target_property(BITPIT_RELEASE_POSTFIX ${BITPIT_LIBRARY} RELEASE_POSTFIX)
get_target_property(BITPIT_INTERFACE_COMPILE_DEFINITIONS ${BITPIT_LIBRARY} INTERFACE_COMPILE_DEFINITIONS)

set(BITPIT_ENABLED_MODULE_LIST "")
foreach(MODULE_NAME IN LISTS BITPIT_MODULE_LIST)
    isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)
    if(MODULE_ENABLED)
        list(APPEND BITPIT_ENABLED_MODULE_LIST ${MODULE_NAME})
    endif()
endforeach()

configure_file("BITPITConfig.cmake.in"  "${CMAKE_CURRENT_BINARY_DIR}/BITPITConfig.cmake" @ONLY)

unset(BITPIT_ENABLED_MODULE_LIST)
unset(BITPIT_INTERFACE_COMPILE_DEFINITIONS)
unset(BITPIT_RELEASE_POSTFIX)
unset(BITPIT_DEBUG_POSTFIX)
unset(BITPIT_INSTALL_PREFIX_CODE)

#------------------------------------------------------------------------------------#
# Generate BITPITConfigVersion.cmake.in
#------------------------------------------------------------------------------------#

configure_file("BITPITConfigVersion.cmake.in"  "${CMAKE_CURRENT_BINARY_DIR}/BITPITConfigVersion.cmake" @ONLY)

#------------------------------------------------------------------------------------#
# CMake targets
#------------------------------------------------------------------------------------#
add_custom_target(clean-cmake COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/BITPITConfig.cmake" DESTINATION "${BITPIT_INSTALL_CMAKEDIR}")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/BITPITConfigVersion.cmake" DESTINATION "${BITPIT_INSTALL_CMAKEDIR}")

install(FILES "FindBITPIT.cmake" DESTINATION "${BITPIT_INSTALL_CMAKEDIR}")
install(FILES "UseBITPIT.cmake" DESTINATION "${BITPIT_INSTALL_CMAKEDIR}")
install(FILES "LibFindMacros.cmake" DESTINATION "${BITPIT_INSTALL_CMAKEDIR}")
