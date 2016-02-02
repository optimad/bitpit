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

#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

# Add a target to generate the examples
option(BUILD_EXAMPLES "Create and install the examples" OFF)

include_directories("${PROJECT_SOURCE_DIR}/src/common")
include_directories("${PROJECT_SOURCE_DIR}/src/operators")
include_directories("${PROJECT_SOURCE_DIR}/src/PABLO")

if(BUILD_EXAMPLES)
	# List of examples
	set(EXAMPLE_LIST "")
	list(APPEND EXAMPLE_LIST "operators")
	list(APPEND EXAMPLE_LIST "PABLO_bubbles_2D")
	list(APPEND EXAMPLE_LIST "PABLO_bubbles_3D")
	list(APPEND EXAMPLE_LIST "PABLO_example_00001")
	list(APPEND EXAMPLE_LIST "PABLO_example_00002")
	list(APPEND EXAMPLE_LIST "PABLO_example_00003")
	list(APPEND EXAMPLE_LIST "PABLO_example_00004")
	list(APPEND EXAMPLE_LIST "PABLO_example_00005")
	list(APPEND EXAMPLE_LIST "PABLO_example_00006")
	list(APPEND EXAMPLE_LIST "PABLO_example_00007")
	list(APPEND EXAMPLE_LIST "PABLO_example_00008")
	list(APPEND EXAMPLE_LIST "PABLO_example_00009")
	list(APPEND EXAMPLE_LIST "PABLO_example_00010")
	list(APPEND EXAMPLE_LIST "PABLO_example_00011")
	list(APPEND EXAMPLE_LIST "PABLO_example_00012")

	#Rules to build the examples
	foreach(EXAMPLE_NAME IN LISTS EXAMPLE_LIST)
		set(EXAMPLE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/${EXAMPLE_NAME}.cpp")
		list(APPEND EXAMPLES_SOURCES ${EXAMPLE_SOURCES})

		add_executable(${EXAMPLE_NAME} "${EXAMPLE_SOURCES}")
		target_link_libraries(${EXAMPLE_NAME} ${BITPIT_LIBRARY})
		install(TARGETS ${EXAMPLE_NAME} DESTINATION doc/${CMAKE_PROJECT_NAME}-${BITPIT_VERSION}/examples)
	endforeach()

	add_custom_target(examples DEPENDS ${EXAMPLE_LIST})
	add_custom_target(clean-examples COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

	set (EXAMPLES_SOURCES ${EXAMPLES_SOURCES} PARENT_SCOPE)
endif()
