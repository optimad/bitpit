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
list(APPEND TESTS "test_levelset_00001")
list(APPEND TESTS "test_levelset_00002")
list(APPEND TESTS "test_levelset_00003")
list(APPEND TESTS "test_levelset_00004")
list(APPEND TESTS "test_levelset_00005")
list(APPEND TESTS "test_levelset_00006")
list(APPEND TESTS "test_levelset_00007")
list(APPEND TESTS "test_levelset_00008")
list(APPEND TESTS "test_levelset_00009")
list(APPEND TESTS "test_levelset_00010")
if (BITPIT_ENABLE_MPI)
    list(APPEND TESTS "test_levelset_parallel_00001:3")
    list(APPEND TESTS "test_levelset_parallel_00002:3")
    list(APPEND TESTS "test_levelset_parallel_00003:3")
endif()

# Test extra modules
set(TEST_EXTRA_MODULES "")

# Test extra libraries
set(TEST_EXTRA_LIBRARIES "")

# Add tests
addModuleIntegrationTests(${MODULE_NAME} "${TESTS}" "${TEST_EXTRA_MODULES}" "${TEST_EXTRA_LIBRARIES}")
unset(TESTS)

# Data files
add_custom_command(
    TARGET "integration_test_levelset_00001" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/naca0012_coarse.dgf" "${CMAKE_CURRENT_BINARY_DIR}/data/naca0012_coarse.dgf"
)

add_custom_command(
    TARGET "integration_test_levelset_00002" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/cube.stl" "${CMAKE_CURRENT_BINARY_DIR}/data/cube.stl"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/sphere.stl" "${CMAKE_CURRENT_BINARY_DIR}/data/sphere.stl"
)

add_custom_command(
    TARGET "integration_test_levelset_00003" PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/naca0012.dgf" "${CMAKE_CURRENT_BINARY_DIR}/data/naca0012.dgf"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/square.dgf" "${CMAKE_CURRENT_BINARY_DIR}/data/square.dgf"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/rectangle.dgf" "${CMAKE_CURRENT_BINARY_DIR}/data/rectangle.dgf"
)

if (BITPIT_ENABLE_MPI)
    add_custom_command(
        TARGET "integration_test_levelset_parallel_00001" PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/cube.stl" "${CMAKE_CURRENT_BINARY_DIR}/data/cube.stl"
    )
    add_custom_command(
        TARGET "integration_test_levelset_parallel_00002" PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/naca0012.dgf" "${CMAKE_CURRENT_BINARY_DIR}/data/naca0012.dgf"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/square.dgf" "${CMAKE_CURRENT_BINARY_DIR}/data/square.dgf"
    )
endif ()
