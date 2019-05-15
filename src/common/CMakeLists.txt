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
#----------------------------------------------------------------------------

# CMake settings
cmake_minimum_required(VERSION 2.8)
initialize_lto_property()

# Name of the current module
get_filename_component(MODULE_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)

# Reset the variable
set(HEADERS "")

# Generate the bitpit header file
set(INCLUDE_LIST "#include \"bitpit_common.hpp\"\n")
foreach(MODULE_NAME IN LISTS BITPIT_MODULE_LIST)
	isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)
	if (MODULE_ENABLED)
		set(INCLUDE_LIST "${INCLUDE_LIST}#include \"bitpit_${MODULE_NAME}.hpp\"\n")
	endif ()
endforeach ()

set(BITPIT_HEADER_TEMPLATE "bitpit_template.hpp.in")
set(BITPIT_HEADER "bitpit.hpp")
CONFIGURE_FILE("${BITPIT_HEADER_TEMPLATE}" "${BITPIT_HEADER}")
list(APPEND HEADERS "${CMAKE_CURRENT_BINARY_DIR}/${BITPIT_HEADER}")

unset(BITPIT_HEADER_TEMPLATE)
unset(BITPIT_HEADER)
unset(INCLUDE_LIST)

# Add common files to targets
addModuleIncludeDirectories(${MODULE_NAME})

file(GLOB SOURCE_FILES "*.cpp")
set(COMMON_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Common sources" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
list(APPEND HEADERS "${HEADER_FILES}")
unset(HEADER_FILES)
set(COMMON_HEADERS "${HEADERS}" CACHE INTERNAL "Common headers" FORCE)

if (NOT "${COMMON_SOURCES}" STREQUAL "")
    set(COMMON_TARGET_OBJECT "COMMON_TARGET_OBJECT")
    add_library(${COMMON_TARGET_OBJECT} OBJECT ${COMMON_SOURCES})
    set_lto_property(${COMMON_TARGET_OBJECT})
endif ()

install(FILES ${COMMON_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
