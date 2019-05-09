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

# Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

# Name of the current module
get_filename_component(MODULE_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)

# List of tests
set(TESTS "")
list(APPEND TESTS "test_levelset_00001")
list(APPEND TESTS "test_levelset_00002")
list(APPEND TESTS "test_levelset_00003")
if (ENABLE_MPI)
	list(APPEND TESTS "test_levelset_parallel_00001:3")
	list(APPEND TESTS "test_levelset_parallel_00002:3")
	list(APPEND TESTS "test_levelset_parallel_00003:3")
endif()

# Test extra libraries
set(TEST_EXTRA_LIBRARIES "")

# Add tests
addModuleTests(${MODULE_NAME} "${TESTS}" "${TEST_EXTRA_LIBRARIES}")
unset(TESTS)

# Data files
add_custom_command(
    TARGET "test_levelset_00002" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/cube.stl" "${CMAKE_CURRENT_BINARY_DIR}/data/cube.stl"
)

add_custom_command(
    TARGET "test_levelset_00003" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/naca0012.dgf" "${CMAKE_CURRENT_BINARY_DIR}/data/naca0012.dgf"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/square.dgf" "${CMAKE_CURRENT_BINARY_DIR}/data/square.dgf"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/rectangle.dgf" "${CMAKE_CURRENT_BINARY_DIR}/data/rectangle.dgf"
)

if (ENABLE_MPI)
	add_custom_command(
		TARGET "test_levelset_parallel_00001" PRE_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/cube.stl" "${CMAKE_CURRENT_BINARY_DIR}/data/cube.stl"
	)
	add_custom_command(
		TARGET "test_levelset_parallel_00002" PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/naca0012.dgf" "${CMAKE_CURRENT_BINARY_DIR}/data/naca0012.dgf"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/square.dgf" "${CMAKE_CURRENT_BINARY_DIR}/data/square.dgf"
	)
endif ()
