#---------------------------------------------------------------------------
#
#  bitpit
#
#  Copyright (C) 2015-2019 OPTIMAD engineering Srl
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
# CMake settings
#------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 2.8)

#------------------------------------------------------------------------------------#
# Subdirectories
#------------------------------------------------------------------------------------#

# Modules
foreach(MODULE_NAME IN LISTS BITPIT_MODULE_LIST)
    isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)
    if (MODULE_ENABLED)
        if (IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE_NAME}")
            add_subdirectory(${MODULE_NAME})
        endif ()
    endif ()
endforeach()

#------------------------------------------------------------------------------------#
# Targets
#------------------------------------------------------------------------------------#
add_custom_target(integration-tests DEPENDS ${TEST_TARGETS})
add_custom_target(clean-integration-tests COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
