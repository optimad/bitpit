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

#------------------------------------------------------------------------------------#
# CMake settings
#------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.10)

#------------------------------------------------------------------------------------#
# Variables
#------------------------------------------------------------------------------------#
set(BITPIT_SOURCE_TEST_DIR ${CMAKE_CURRENT_LIST_DIR})

#------------------------------------------------------------------------------------#
# Functions
#------------------------------------------------------------------------------------#

# Get the name of the test
function(getTestName TEST_ENTRY TEST_NAME)
    string(REGEX REPLACE ":[^:]*$" "" TEST_ENTRY_STRIPPED ${TEST_ENTRY})
    set(${TEST_NAME} "${TEST_ENTRY_STRIPPED}" PARENT_SCOPE)
endfunction()

# Check whether the test is parallel or not
function(isTestParallel TEST_ENTRY PARALLEL_FLAG)
    if (TEST_ENTRY MATCHES "parallel")
        set(${PARALLEL_FLAG} "1" PARENT_SCOPE)
    else ()
        set(${PARALLEL_FLAG} "0" PARENT_SCOPE)
    endif ()
endfunction()

# Get the number of processes a test has to run on
function(getTestNumProcs TEST_ENTRY N_PROCS)
    isTestParallel("${TEST_ENTRY}" IS_TEST_PARALLEL)
    if (NOT IS_TEST_PARALLEL)
        set(${N_PROCS} "1" PARENT_SCOPE)
        return()
    endif ()

    string(REGEX MATCH "[^:]*$" TEST_ENTRY_MATCH ${TEST_ENTRY})
    if ("${TEST_ENTRY_MATCH}" STREQUAL "${TEST_ENTRY}")
        set(BITPIT_TEST_DEFAULT_NUMPROCS 2 CACHE STRING "Is the default number of precesses the test will run on")
        mark_as_advanced(BITPIT_TEST_DEFAULT_NUMPROCS)
        set(${N_PROCS} "${BITPIT_TEST_DEFAULT_NUMPROCS}" PARENT_SCOPE)
        return()
    endif ()

    set(${N_PROCS} "${TEST_ENTRY_MATCH}" PARENT_SCOPE)
endfunction()

# Add an integration test for a module
function(addModuleIntegrationTests MODULE_NAME TEST_ENTRIES TEST_LIBRARIES)
    addModuleTests(${MODULE_NAME} "integration" "${TEST_ENTRIES}" "${TEST_LIBRARIES}")
endfunction()

# Add a unit test for a module
function(addModuleUnitTests MODULE_NAME TEST_ENTRIES TEST_LIBRARIES)
    addModuleTests(${MODULE_NAME} "unit" "${TEST_ENTRIES}" "${TEST_LIBRARIES}")
endfunction()

# Add a test for a module
function(addModuleTests MODULE_NAME TEST_TYPE TEST_ENTRIES TEST_LIBRARIES)
    # Add single tests
    set(TARGETS "")
    foreach (TEST_ENTRY IN LISTS TEST_ENTRIES)
        isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)
        if (NOT MODULE_ENABLED)
            return ()
        endif ()

        isTestParallel("${TEST_ENTRY}" IS_TEST_PARALLEL)
        if (IS_TEST_PARALLEL AND NOT BITPIT_ENABLE_MPI)
            return ()
        endif ()

        getTestName("${TEST_ENTRY}" TEST_NAME)
        if (${IS_TEST_PARALLEL})
            getTestNumProcs("${TEST_ENTRY}" N_PROCS)
            addModuleParallelTest("${TEST_NAME}" "${TEST_TYPE}" "${TEST_LIBRARIES}" ${N_PROCS} TEST_TARGET)
        else ()
            addModuleSerialTest("${TEST_NAME}" "${TEST_TYPE}" "${TEST_LIBRARIES}" TEST_TARGET)
        endif ()
        list(APPEND TARGETS "${TEST_TARGET}")
    endforeach ()

    list(LENGTH TARGETS TESTS_COUNT)
    if (${TESTS_COUNT} EQUAL 0)
        return()
    endif ()

    UNSET(MODULE_HAS_TESTS)

    # Test targets for the module
    string(TOUPPER ${MODULE_NAME} UPPERCASE_MODULE_NAME)
    set(${UPPERCASE_MODULE_NAME}_TEST_TARGETS "${TARGETS}" CACHE INTERNAL "List of tests targets for the ${MODULE_NAME} module" FORCE)

    # Add the targets to the global list of targets
    set(TMP_TEST_TARGETS "${TEST_TARGETS}")
    foreach (TARGET IN LISTS TARGETS)
        list(APPEND TMP_TEST_TARGETS "${TARGET}")
    endforeach ()
    set(TEST_TARGETS "${TMP_TEST_TARGETS}" CACHE INTERNAL "List of tests targets" FORCE)

    # Add include directories required by the module
    addModuleIncludeDirectories(${MODULE_NAME})

    # Add include directories required by the test framework
    include_directories("${BITPIT_SOURCE_TEST_DIR}" "${BITPIT_SOURCE_TEST_DIR}/${TEST_TYPE}_tests")

    # Add rules for the module
    add_custom_target(${TEST_TYPE}-tests-${MODULE_NAME} DEPENDS ${TARGETS})
    add_custom_target(clean-${TEST_TYPE}-tests-${MODULE_NAME} COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

    UNSET(TEST_NAMES)
endfunction()

# Add a serial test
function(addModuleSerialTest TEST_NAME TEST_TYPE TEST_LIBRARIES PARENT_TEST_TARGET_NAME)
    addTest("${TEST_NAME}" "${TEST_TYPE}" "${TEST_LIBRARIES}" "${CMAKE_CURRENT_BINARY_DIR}" 1 TEST_TARGET_NAME)
    set(${PARENT_TEST_TARGET_NAME} "${TEST_TARGET_NAME}" PARENT_SCOPE)
endfunction()

# Add a parallel test
function(addModuleParallelTest TEST_NAME TEST_TYPE TEST_LIBRARIES N_PROCS PARENT_TEST_TARGET_NAME)
    # Parallel tests are enable only in MPI builds
    if (NOT BITPIT_ENABLE_MPI)
        return()
    endif ()

    # Add the test
    addTest("${TEST_NAME}" "${TEST_TYPE}" "${TEST_LIBRARIES}" "${CMAKE_CURRENT_BINARY_DIR}" ${N_PROCS} TEST_TARGET_NAME)
    set(${PARENT_TEST_TARGET_NAME} "${TEST_TARGET_NAME}" PARENT_SCOPE)
endfunction()

# Add a test
function(addTest TEST_NAME TEST_TYPE TEST_LIBRARIES WORKING_DIRECTORY N_PROCS PARENT_TEST_TARGET_NAME)

    # Parallel tests are enable only in MPI builds
    if (NOT BITPIT_ENABLE_MPI AND ${N_PROCS} GREATER 1)
        return()
    endif ()

    # Get name of test target
    set(TEST_TARGET_NAME "${TEST_TYPE}_${TEST_NAME}")
    set(${PARENT_TEST_TARGET_NAME} "${TEST_TARGET_NAME}" PARENT_SCOPE)

    # Test command
    if (BITPIT_ENABLE_MPI)
        set(TEST_EXEC ${MPIEXEC})
        set(TEST_ARGS ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${N_PROCS} ${MPIEXEC_POSTFLAGS} "$<TARGET_FILE:${TEST_TARGET_NAME}>")
    else()
        set(TEST_EXEC "$<TARGET_FILE:${TEST_TARGET_NAME}>")
        set(TEST_ARGS "")
    endif()


    # Add test target
    add_executable(${TEST_TARGET_NAME} "${TEST_NAME}.cpp")
    target_link_libraries(${TEST_TARGET_NAME} ${BITPIT_LIBRARY})
    target_link_libraries(${TEST_TARGET_NAME} ${BITPIT_EXTERNAL_LIBRARIES})
    target_link_libraries(${TEST_TARGET_NAME} ${TEST_LIBRARIES})
    if (CMAKE_BUILD_TYPE STREQUAL "Coverage")
        target_link_libraries(${TEST_TARGET_NAME} "gcov")
    endif()

    # Configure source formatting
    set(TEST_CLANG_FORMAT_TARGET "clang-format-${TEST_TARGET_NAME}")
    ClangFormatAddTarget("${TEST_CLANG_FORMAT_TARGET}")
    ClangFormatTargetSources("${TEST_CLANG_FORMAT_TARGET}" "${TEST_NAME}.cpp")
    add_dependencies(${CLANG_FORMAT_TARGET} "${TEST_CLANG_FORMAT_TARGET}")

    # Add test
    add_test(NAME "${TEST_TARGET_NAME}" COMMAND ${TEST_EXEC} ${TEST_ARGS} WORKING_DIRECTORY "${WORKING_DIRECTORY}")

    # Set test variables
    if (${N_PROCS} GREATER 1)
        set_tests_properties(${TEST_TARGET_NAME} PROPERTIES ENVIRONMENT OMPI_MCA_rmaps_base_oversubscribe=1)
    endif()

endfunction()

#------------------------------------------------------------------------------------#
# Subdirectories
#------------------------------------------------------------------------------------#
set(TEST_TARGETS "" CACHE INTERNAL "List of tests targets" FORCE)

if(BITPIT_ENABLE_INTEGRATION_TESTS)
    add_subdirectory(integration_tests)
endif()

if(BITPIT_ENABLE_UNIT_TESTS)
    add_subdirectory(unit_tests)
endif()

#------------------------------------------------------------------------------------#
# Targets
#------------------------------------------------------------------------------------#
add_custom_target(tests DEPENDS ${TEST_TARGETS})
add_custom_target(clean-tests COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(check DEPENDS tests COMMAND ${CMAKE_MAKE_PROGRAM} test)
