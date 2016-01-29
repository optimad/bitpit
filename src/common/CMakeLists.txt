#---------------------------------------------------------------------------
#
#  bitpit
#
#  Copyright (C) 2015-2016 OPTIMAD engineering Srl
#
#  -------------------------------------------------------------------------
#  License
#  This file is part of bitbit.
#
#  bitpit is free software: you can redistribute it and/or modify it
#  under the terms of the GNU Lesser General Public License v3 (LGPL)
#  as published by by the Free Software Foundation.
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

# Reset the variable
set(HEADERS "")

# Generate the reference header files
set(HEADER_TEMPLATE "bitpit_template.hpp.in")
foreach (GROUP_NAME IN LISTS BITPIT_GROUP_LIST)
	set(REFERENCE_HEADER "bitpit_${GROUP_NAME}.hpp")
	string(TOUPPER ${GROUP_NAME} UPPER_GROUP_NAME)

	set(INCLUDE_LIST "#include <bitpit_Base_common.hpp>\n")
	foreach(MODULE_NAME IN LISTS BITPIT_${UPPER_GROUP_NAME}_MODULE_LIST)
		isModuledEnabled(${MODULE_NAME} MODULE_ENABLED)
		if (MODULE_ENABLED)
			set(INCLUDE_LIST "${INCLUDE_LIST}#include <bitpit_Base_${MODULE_NAME}.hpp>\n")
		endif ()
	endforeach ()

	CONFIGURE_FILE("${HEADER_TEMPLATE}" "${REFERENCE_HEADER}")
	list(APPEND HEADERS "${CMAKE_CURRENT_BINARY_DIR}/${REFERENCE_HEADER}")

	unset(INCLUDE_LIST)
	unset(REFERENCE_HEADER)
endforeach ()

# Add common files to targets
file(GLOB SOURCE_FILES "*.cpp")
set(COMMON_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Common sources" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
list(APPEND HEADERS "${HEADER_FILES}")
unset(HEADER_FILES)
set(COMMON_HEADERS "${HEADERS}" CACHE INTERNAL "Common headers" FORCE)

if (NOT "${COMMON_SOURCES}" STREQUAL "")
	add_library(COMMON_TARGET_OBJECT OBJECT ${COMMON_SOURCES})
endif ()

install(FILES ${COMMON_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
