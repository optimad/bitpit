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

# Name of the current module
get_filename_component(MODULE_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)

# List of tests
set(TESTS "")
list(APPEND TESTS "test_volcartesian_00001")
list(APPEND TESTS "test_volcartesian_00002")
list(APPEND TESTS "test_volcartesian_00003")

set(VOLCARTESIAN_TEST_ENTRIES "${TESTS}" CACHE INTERNAL "List of tests for the volcartesian module" FORCE)

# Rules to build the tests
addModuleIncludeDirectories(${MODULE_NAME})

set(TEST_TARGETS "")
foreach(TEST_ENTRY IN LISTS VOLCARTESIAN_TEST_ENTRIES)
	getTestName("${TEST_ENTRY}" TEST_NAME)
	set(TEST_SOURCES "${TEST_NAME}.cpp")
	list(APPEND TEST_TARGETS "${TEST_NAME}")

	add_executable(${TEST_NAME} "${TEST_SOURCES}")
	target_link_libraries(${TEST_NAME} ${BITPIT_LIBRARY})
	target_link_libraries(${TEST_NAME} ${BITPIT_EXTERNAL_LIBRARIES})
endforeach()
set(VOLCARTESIAN_TEST_TARGETS "${TEST_TARGETS}" CACHE INTERNAL "List of test targets for the volcartesian module" FORCE)

add_custom_target(tests-volcartesian DEPENDS ${TEST_TARGETS})
add_custom_target(clean-tests-volcartesian COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

unset(TEST_TARGETS)