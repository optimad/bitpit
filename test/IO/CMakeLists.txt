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
list(APPEND TESTS "test_IO_00001")
list(APPEND TESTS "test_IO_00002")

set(IO_TEST_ENTRIES "${TESTS}" CACHE INTERNAL "List of tests fo the IO module" FORCE)

# External dependencies

# Rules to build the tests
include_directories("${PROJECT_SOURCE_DIR}/src/common")
include_directories("${PROJECT_SOURCE_DIR}/src/operators")
include_directories("${PROJECT_SOURCE_DIR}/src/containers")
include_directories("${PROJECT_SOURCE_DIR}/src/IO")

set(TEST_TARGETS "")
foreach(TEST_ENTRY IN LISTS IO_TEST_ENTRIES)
	getTestName("${TEST_ENTRY}" TEST_NAME)
	set(TEST_SOURCES "${TEST_NAME}.cpp")
	list(APPEND TEST_TARGETS "${TEST_NAME}")

    add_executable(${TEST_NAME} "${TEST_SOURCES}")
    target_link_libraries(${TEST_NAME} ${BITPIT_LIBRARY})
	target_link_libraries(${TEST_NAME} ${BITPIT_EXTERNAL_DEPENDENCIES})
endforeach()
set(IO_TEST_TARGETS "${TEST_TARGETS}" CACHE INTERNAL "List of test targets for the IO module" FORCE)

add_custom_target(tests-IO DEPENDS ${TEST_TARGETS})
add_custom_target(clean-tests-IO COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

unset(TEST_TARGETS)

# Data files
add_custom_command(
    TARGET "test_IO_00001" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/selection.vtu" "${CMAKE_CURRENT_BINARY_DIR}/data/selection.vtu"
    )
