#---------------------------------------------------------------------------
#
#  bitpit
#
#  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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

# Name of the current module
get_filename_component(MODULE_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)

# List of tests
set(TESTS "")
list(APPEND TESTS "test_surfunstructured_00001")
list(APPEND TESTS "test_surfunstructured_00002")
list(APPEND TESTS "test_surfunstructured_00003")
list(APPEND TESTS "test_surfunstructured_00004")
list(APPEND TESTS "test_surfunstructured_00005")
list(APPEND TESTS "test_surfunstructured_00006")
list(APPEND TESTS "test_surfunstructured_00007")
list(APPEND TESTS "test_surfunstructured_00008")
list(APPEND TESTS "test_surfunstructured_00009")
if (BITPIT_ENABLE_MPI)
    list(APPEND TESTS "test_surfunstructured_parallel_00001:4")
    list(APPEND TESTS "test_surfunstructured_parallel_00002:2")
    list(APPEND TESTS "test_surfunstructured_parallel_00003:4")
    list(APPEND TESTS "test_surfunstructured_parallel_00004:4")
    list(APPEND TESTS "test_surfunstructured_parallel_00005:3")
    list(APPEND TESTS "test_surfunstructured_parallel_00006:2")
    list(APPEND TESTS "test_surfunstructured_parallel_00007:4")
    list(APPEND TESTS "test_surfunstructured_parallel_00008:3")
endif ()

# Test extra modules
set(TEST_EXTRA_MODULES "")

# Test extra libraries
set(TEST_EXTRA_LIBRARIES "")

# Add tests
addModuleIntegrationTests(${MODULE_NAME} "${TESTS}" "${TEST_EXTRA_MODULES}" "${TEST_EXTRA_LIBRARIES}")
unset(TESTS)

# Data files
add_custom_command(
    TARGET "integration_test_surfunstructured_00003" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/buddha.stl" "${CMAKE_CURRENT_BINARY_DIR}/data/buddha.stl"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/cube.stl" "${CMAKE_CURRENT_BINARY_DIR}/data/cube.stl"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/NACA0012.dgf" "${CMAKE_CURRENT_BINARY_DIR}/data/NACA0012.dgf"
    )

add_custom_command(
    TARGET "integration_test_surfunstructured_00004" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/buddha.stl" "${CMAKE_CURRENT_BINARY_DIR}/data/buddha.stl"
    )

add_custom_command(
    TARGET "integration_test_surfunstructured_00005" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/ahmed.stl" "${CMAKE_CURRENT_BINARY_DIR}/data/ahmed.stl"
    )

add_custom_command(
    TARGET "integration_test_surfunstructured_00009" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/buddha.stl" "${CMAKE_CURRENT_BINARY_DIR}/data/buddha.stl"
    )

if (BITPIT_ENABLE_MPI)
    add_custom_command(
        TARGET "integration_test_surfunstructured_parallel_00002" PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/buddha.stl" "${CMAKE_CURRENT_BINARY_DIR}/data/buddha.stl"
    )

    add_custom_command(
        TARGET "integration_test_surfunstructured_parallel_00003" PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/buddha.stl" "${CMAKE_CURRENT_BINARY_DIR}/data/buddha.stl"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/disc.stl" "${CMAKE_CURRENT_BINARY_DIR}/data/disc.stl"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/circle.dgf" "${CMAKE_CURRENT_BINARY_DIR}/data/circle.dgf"
    )

    add_custom_command(
        TARGET "integration_test_surfunstructured_parallel_00005" PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/sphere.stl" "${CMAKE_CURRENT_BINARY_DIR}/data/sphere.stl"
    )
endif ()
