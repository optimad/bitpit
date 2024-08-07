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
list(APPEND TESTS "test_LA_00001")
list(APPEND TESTS "test_LA_00002")
list(APPEND TESTS "test_LA_00003")
list(APPEND TESTS "test_LA_00004")
list(APPEND TESTS "test_LA_00005")
list(APPEND TESTS "test_LA_00006")
list(APPEND TESTS "test_LA_00007")
if (BITPIT_ENABLE_MPI)
    list(APPEND TESTS "test_LA_parallel_00001")
    list(APPEND TESTS "test_LA_parallel_00002")
    list(APPEND TESTS "test_LA_parallel_00003")
    list(APPEND TESTS "test_LA_parallel_00004")
    list(APPEND TESTS "test_LA_parallel_00005")
endif()

# Test extra modules
set(TEST_EXTRA_MODULES "")

# Test extra libraries
set(TEST_EXTRA_LIBRARIES "")

# Add tests
addModuleIntegrationTests(${MODULE_NAME} "${TESTS}" "${TEST_EXTRA_MODULES}" "${TEST_EXTRA_LIBRARIES}")
unset(TESTS)
