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

# Name of the current module
get_filename_component(MODULE_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)

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
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${BITPIT_HEADER}" DESTINATION include/${PROJECT_NAME})

unset(BITPIT_HEADER_TEMPLATE)
unset(BITPIT_HEADER)
unset(INCLUDE_LIST)


# Configure module
configureModule(${MODULE_NAME})
