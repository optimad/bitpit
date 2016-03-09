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

# Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

# List of tests
set(TESTS "")
list(APPEND TESTS "test_surftripatch_00001")
list(APPEND TESTS "test_surftripatch_00002")
list(APPEND TESTS "test_surftripatch_00003")
list(APPEND TESTS "test_surftripatch_00004")
if (ENABLE_MPI)
	list(APPEND TESTS "test_surftripatch_parallel_00001:4")
	list(APPEND TESTS "test_surftripatch_parallel_00002:2")
endif ()

set(SURFTRIPATCH_TEST_ENTRIES "${TESTS}" CACHE INTERNAL "List of tests for the surftripatch module" FORCE)

# Rules to build the tests
include_directories("${PROJECT_SOURCE_DIR}/src/common")
include_directories("${PROJECT_SOURCE_DIR}/src/operators")
include_directories("${PROJECT_SOURCE_DIR}/src/containers")
include_directories("${PROJECT_SOURCE_DIR}/src/IO")
include_directories("${PROJECT_SOURCE_DIR}/src/patchkernel")
include_directories("${PROJECT_SOURCE_DIR}/src/surftripatch")

set(TEST_TARGETS "")
foreach(TEST_ENTRY IN LISTS SURFTRIPATCH_TEST_ENTRIES)
	getTestName("${TEST_ENTRY}" TEST_NAME)
	set(TEST_SOURCES "${TEST_NAME}.cpp")
	list(APPEND TEST_TARGETS "${TEST_NAME}")

	add_executable(${TEST_NAME} "${TEST_SOURCES}")
	target_link_libraries(${TEST_NAME} ${BITPIT_LIBRARY})
endforeach()
set(SURFTRIPATCH_TEST_TARGETS "${TEST_TARGETS}" CACHE INTERNAL "List of test targets for the surftripatch module" FORCE)

add_custom_target(tests-surftripatch DEPENDS ${TEST_TARGETS})
add_custom_target(clean-tests-surftripatch COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

unset(TEST_TARGETS)

# Data files
add_custom_command(
    TARGET "test_surftripatch_00003" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/buddha.stl" "${CMAKE_CURRENT_BINARY_DIR}/data/buddha.stl"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/cube.stl" "${CMAKE_CURRENT_BINARY_DIR}/data/cube.stl"
    )

add_custom_command(
    TARGET "test_surftripatch_00004" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/buddha.stl" "${CMAKE_CURRENT_BINARY_DIR}/data/buddha.stl"
    )

if (ENABLE_MPI)
	add_custom_command(
		TARGET "test_surftripatch_parallel_00002" PRE_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/buddha.stl" "${CMAKE_CURRENT_BINARY_DIR}/data/buddha.stl"
	)
endif ()
