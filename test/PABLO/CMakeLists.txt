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

# List serial of tests
set(TESTS_SERIAL "")
list(APPEND TESTS_SERIAL "test_PABLO_00001")
list(APPEND TESTS_SERIAL "test_PABLO_00002")
list(APPEND TESTS_SERIAL "test_PABLO_00003")
list(APPEND TESTS_SERIAL "test_PABLO_00004")

set(PABLO_TEST_LIST_SERIAL "${TESTS_SERIAL}" CACHE INTERNAL "List of serial tests for the PABLO module" FORCE)

# List of parallel tests
if (ENABLE_MPI)
	set(TESTS_PARALLEL "")
	list(APPEND TESTS_PARALLEL "test_PABLO_parallel_00001")

	set(PABLO_TEST_LIST_PARALLEL "${TESTS_PARALLEL}" CACHE INTERNAL "List of parallel tests for the PABLO module" FORCE)
endif()

# Rules to build the tests
include_directories("${PROJECT_SOURCE_DIR}/src/operators")
include_directories("${PROJECT_SOURCE_DIR}/src/containers")
include_directories("${PROJECT_SOURCE_DIR}/src/PABLO")

foreach(TEST_NAME IN LISTS PABLO_TEST_LIST_SERIAL PABLO_TEST_LIST_PARALLEL)
	set(TEST_SOURCES "${TEST_NAME}.cpp")

	add_executable(${TEST_NAME} "${TEST_SOURCES}")
	target_link_libraries(${TEST_NAME} ${BITPIT_LIBRARY})
endforeach()

add_custom_target(tests-PABLO DEPENDS ${PABLO_TEST_LIST_SERIAL})
add_custom_target(clean-tests-PABLO COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
