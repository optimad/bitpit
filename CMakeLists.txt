# ---------------------------------------------------------------------------
#
# bitpit
#
# Copyright (C) 2015-2021 OPTIMAD engineering Srl
#
# -------------------------------------------------------------------------
# License
# This file is part of bitpit.
#
# bitpit is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License v3 (LGPL)
# as published by the Free Software Foundation.
#
# bitpit is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with bitpit. If not, see <http://www.gnu.org/licenses/>.
#
# ---------------------------------------------------------------------------*/

# ------------------------------------------------------------------------------------#
# CMake settings
# ------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.13)

project("bitpit" CXX)

# ------------------------------------------------------------------------------------#
# External utilities
# ------------------------------------------------------------------------------------#
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/external/git/cmake")
include(GetGitRevisionDescription)

include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/external/clang-format/cmake")
include(ClangFormat)

include(GNUInstallDirs)

include(CTest)

# ------------------------------------------------------------------------------------#
# Variables visible to the user
#
# NOTE: to generate Windows lib files, both BUILD_SHARED_LIBS and CMAKE_GNUtoMS flags
# should be enabled.
# ------------------------------------------------------------------------------------#
set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries")

set(BITPIT_ENABLE_DOC OFF CACHE BOOL "If set, the HTML-based API documentation will be created (requires Doxygen)")
set(BITPIT_ENABLE_MPI ON CACHE BOOL "If set, the program is compiled with MPI support")

set(BITPIT_LTO_STRATEGY "Auto" CACHE STRING "Choose the Link Time Optimization (LTO) strategy, options are: Auto (i.e., optimiziation is enabled only in release build and only for some tested configurations) Enabled Disabled.")
set_property(CACHE BITPIT_LTO_STRATEGY PROPERTY STRINGS "Auto" "Enabled" "Disabled")

# ------------------------------------------------------------------------------------#
# Internal variables
# ------------------------------------------------------------------------------------#
set(BITPIT_VERSION "1.9.0-devel" CACHE INTERNAL "Version number")
set(BITPIT_INSTALL_CMAKEDIR "${CMAKE_INSTALL_DATAROOTDIR}/${CMAKE_PROJECT_NAME}/cmake" CACHE INTERNAL "Directory in which cmake modules will be installed")

if(MSVC)
    # export bitpit all symbols definition in a *.def file when create the *.dll dynamic library
    # and its import lib for target linking (dynamic *.lib). Def file will be needed when an external
    # project (library or executable) will need bitpit.
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

# ------------------------------------------------------------------------------------#
# Functions
# ------------------------------------------------------------------------------------#
function(getModuleEnableName MODULE_NAME FLAG_NAME)
    string(TOUPPER ${MODULE_NAME} UPPER_MODULE_NAME)
    set(${FLAG_NAME} "BITPIT_MODULE_${UPPER_MODULE_NAME}" PARENT_SCOPE)
endfunction()

function(isModuleEnabled MODULE_NAME ENABLED)
    getModuleEnableName(${MODULE_NAME} ENABLED_VARIABLE)

    if(DEFINED ${ENABLED_VARIABLE})
        set(${ENABLED} ${${ENABLED_VARIABLE}} PARENT_SCOPE)
    else()
        set(${ENABLED} "FALSE" PARENT_SCOPE)
    endif()
endfunction()

function(areModulesEnabled MODULE_LIST ENABLED)
    foreach(MODULE_NAME IN LISTS MODULE_LIST)
        isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)

        if(NOT MODULE_ENABLED)
            set(${ENABLED} "FALSE" PARENT_SCOPE)
            return()
        endif()
    endforeach()

    set(${ENABLED} "TRUE" PARENT_SCOPE)
endfunction()

function(enableModule MODULE_NAME ENABLED)
    set(EXTRA_ARGUMENTS ${ARGN})
    list(LENGTH EXTRA_ARGUMENTS EXTRA_ARGUMENT_COUNT)

    if(${EXTRA_ARGUMENT_COUNT} GREATER 0)
        list(GET EXTRA_ARGUMENTS 0 FORCED)

        if(FORCED)
            set(FORCE_FLAG "FORCE")
        endif()
    endif()

    if(NOT DEFINED FORCE_FLAG)
        set(FORCE_FLAG "")
    endif()

    getModuleEnableName(${MODULE_NAME} MODULE_ENABLE_FLAG)
    set(${MODULE_ENABLE_FLAG} ${ENABLED} CACHE BOOL "Request building ${MODULE_NAME} module" ${FORCE_FLAG})
    mark_as_advanced(${MODULE_ENABLE_FLAG})
endfunction()

function(hideModule MODULE_NAME)
    getModuleEnableName(${MODULE_NAME} MODULE_ENABLE_FLAG)

    if(DEFINED ${MODULE_ENABLE_FLAG})
        unset(${MODULE_ENABLE_FLAG} CACHE)
    endif()
endfunction()

function(addModuleIncludeDirectories TARGET_NAME MODULE_NAME)
    set(EXTRA_ARGUMENTS ${ARGN})
    list(LENGTH EXTRA_ARGUMENTS EXTRA_ARGUMENT_COUNT)

    if(${EXTRA_ARGUMENT_COUNT} GREATER 0)
        list(GET EXTRA_ARGUMENTS 0 INCLUDE_EXTERNAL)
    else()
        set(INCLUDE_EXTERNAL "TRUE")
    endif()

    # External includes
    if(INCLUDE_EXTERNAL)
        target_include_directories(${TARGET_NAME} PUBLIC "${BITPIT_EXTERNAL_INCLUDE_DIRS}")
    endif()

    # Add dependiecies
    string(TOUPPER ${MODULE_NAME} UPPER_MODULE_NAME)

    foreach(DEPENDENCY_NAME IN LISTS ${UPPER_MODULE_NAME}_DEPS)
        addModuleIncludeDirectories(${TARGET_NAME} ${DEPENDENCY_NAME} "FALSE")
    endforeach()

    unset(UPPER_MODULE_NAME)

    # Add module directory
    target_include_directories(${TARGET_NAME} PRIVATE "${PROJECT_SOURCE_DIR}/src/${MODULE_NAME}")
    target_include_directories(${TARGET_NAME} PRIVATE "${PROJECT_BINARY_DIR}/src/${MODULE_NAME}")
endfunction()

# Define an experimental feature
#
# Experimental features are disabled by default.
function(defineExperimentalFeature FEATURE_NAME FEATURE_DESCRIPTION)
    defineFeature("EXPERIMENTAL_${FEATURE_NAME}" ${FEATURE_DESCRIPTION} 0)
endfunction()

# Define a module experimental feature
#
# Experimental features are disabled by default.
function(defineExperimentalModuleFeature MODULE_NAME FEATURE_NAME FEATURE_DESCRIPTION)
    string(TOUPPER ${MODULE_NAME} UPPERCASE_MODULE_NAME)
    set(MODULE_FEATURE_NAME "${UPPERCASE_MODULE_NAME}_EXPERIMENTAL_${FEATURE_NAME}")
    defineFeature("${MODULE_FEATURE_NAME}" ${FEATURE_DESCRIPTION} 0)
endfunction()

# Define a deprecated feature
#
# Deprecated features are enabled by default.
function(defineDeprecatedFeature FEATURE_NAME FEATURE_DESCRIPTION)
    defineFeature("DEPRECATED_${FEATURE_NAME}" ${FEATURE_DESCRIPTION} 1)
endfunction()

# Define a module deprecated feature
#
# Deprecated features are enabled by default.
function(defineDeprecatedModuleFeature MODULE_NAME FEATURE_NAME FEATURE_DESCRIPTION)
    string(TOUPPER ${MODULE_NAME} UPPERCASE_MODULE_NAME)
    set(MODULE_FEATURE_NAME "${UPPERCASE_MODULE_NAME}_DEPRECATED_${FEATURE_NAME}")
    defineFeature("${MODULE_FEATURE_NAME}" ${FEATURE_DESCRIPTION} 1)
endfunction()

# Define a feature
#
# For each feature a corresponding preprocessor macro will be defined. Given
# the feature FEATURE the corresponding preprocessor macro will be called
# BITPIT_FEATURE. Ye value of the macro will be set to 1 if the feature is
# enable or it will be set to 0 if the feature is disabled.
function(defineFeature FEATURE_NAME FEATURE_DESCRIPTION DEFAULT_STATUS)
    string(TOUPPER ${FEATURE_NAME} UPPER_FEATURE_NAME)
    set(BITPIT_FEATURE_NAME BITPIT_${UPPER_FEATURE_NAME})
    set(${BITPIT_FEATURE_NAME} ${DEFAULT_STATUS} CACHE BOOL ${FEATURE_DESCRIPTION})
    mark_as_advanced(${BITPIT_FEATURE_NAME})

    set(FEATURE_DEFINITION "${BITPIT_FEATURE_NAME}")

    if(NOT DEFINED ${BITPIT_FEATURE_NAME})
        set(FEATURE_DEFINITION "${FEATURE_DEFINITION}=0")
    elseif(${BITPIT_FEATURE_NAME})
        set(FEATURE_DEFINITION "${FEATURE_DEFINITION}=1")
    else()
        set(FEATURE_DEFINITION "${FEATURE_DEFINITION}=0")
    endif()

    target_compile_definitions(${BITPIT_LIBRARY} PUBLIC "${FEATURE_DEFINITION}")
endfunction()

# Add coverage analysis
function(addCoverageAnalysis GENERATE_GCOVR_REPORTS GENERATE_LCOV_REPORTS EXCLUDES)
    # Collect excludes
    set(OVERALL_EXCLUDES "/usr/include/*;/usr/lib/*")

    foreach(EXCLUDE ${EXCLUDES} ${COVERAGE_EXCLUDES})
        list(APPEND OVERALL_EXCLUDES "${EXCLUDE}")
    endforeach()

    list(REMOVE_DUPLICATES OVERALL_EXCLUDES)

    # Coverage output directory
    set(COVERAGE_OUTPUT_DIR "${PROJECT_BINARY_DIR}/coverage")

    # Initialize coverage analysis
    #
    # Before initializing coverage analyses needs we need to build the tests.
    add_custom_target(coverage-initialize

        # Remove existing code coverage data
        COMMAND find ${PROJECT_BINARY_DIR} -name "*.gcda" -exec rm {} \;

        VERBATIM
        DEPENDS tests
        COMMENT "Initializeing coverage analysis."
    )

    # Tests are used to gather coverage data
    add_custom_target(coverage
        COMMAND ${CMAKE_COMMAND} -E make_directory ${COVERAGE_OUTPUT_DIR}
        COMMAND ${CMAKE_CTEST_COMMAND}

        BYPRODUCTS ${COVERAGE_OUTPUT_DIR}
        DEPENDS coverage-initialize
        COMMENT "Running tests to gather coverage data."
    )

    # Generate gcovr reports
    if(GENERATE_GCOVR_REPORTS)
        # Find needed executables
        find_program(GCOVR_PATH gcovr REQUIRED)

        # Combine excludes to several -e arguments
        set(GCOVR_EXCLUDE_ARGS "")

        foreach(EXCLUDE ${OVERALL_EXCLUDES})
            list(APPEND GCOVR_EXCLUDE_ARGS "-e")
            list(APPEND GCOVR_EXCLUDE_ARGS "${EXCLUDE}")
        endforeach()

        # Set output directory
        set(GCOVR_OUTPUT_DIR "${COVERAGE_OUTPUT_DIR}/gcovr")

        # Initialize reports
        add_custom_command(TARGET coverage-initialize POST_BUILD

            # Create output directory
            COMMAND ${CMAKE_COMMAND} -E make_directory ${GCOVR_OUTPUT_DIR}

            BYPRODUCTS ${GCOVR_OUTPUT_DIR}
            VERBATIM
            COMMENT "Initializing gcovr reports"
        )

        # Generate reports
        add_custom_command(TARGET coverage POST_BUILD
            COMMAND ${GCOVR_PATH} --html --html-details
            -r ${PROJECT_SOURCE_DIR} ${GCOVR_EXCLUDE_ARGS}
            --object-directory=${PROJECT_BINARY_DIR}
            -o ${GCOVR_OUTPUT_DIR}/index.html

            COMMAND ${GCOVR_PATH} --xml
            -r ${PROJECT_SOURCE_DIR} ${GCOVR_EXCLUDE_ARGS}
            --object-directory=${PROJECT_BINARY_DIR}
            -o ${GCOVR_OUTPUT_DIR}/coverage.xml

            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
            VERBATIM
            COMMENT "Running gcovr to produce code coverage reports."
        )

        # Show path of the reports
        add_custom_command(TARGET coverage POST_BUILD
            COMMAND ;
            COMMENT "gcovr HTML coverage report saved in \"${GCOVR_OUTPUT_DIR}/index.html.html\"."
        )

        add_custom_command(TARGET coverage POST_BUILD
            COMMAND ;
            COMMENT "gcovr XML coverage report saved in \"${GCOVR_OUTPUT_DIR}/coverage.xml\"."
        )
    endif()

    # Generate lcov reports
    if(GENERATE_LCOV_REPORTS)
        # Set output directory
        set(LCOV_OUTPUT_DIR "${COVERAGE_OUTPUT_DIR}/lcov")
        set(LCOV_CASE_NAME "coverage")

        # Find needed executables
        find_program(GCOV_PATH gcov REQUIRED)
        find_program(LCOV_PATH NAMES lcov lcov.bat lcov.exe lcov.perl REQUIRED)
        find_program(GENHTML_PATH NAMES genhtml genhtml.perl genhtml.bat REQUIRED)
        find_program(CPPFILT_PATH NAMES c++filt)

        # lcov arguments
        set(LCOV_ARGS "")

        # genhtml arguments
        set(GENHTML_ARGS "")

        if(CPPFILT_PATH)
            set(GENHTML_EXTRA_ARGS ${GENHTML_ARGS} "--demangle-cpp")
        endif()

        # Initialize lcov reports
        add_custom_command(TARGET coverage-initialize POST_BUILD

            # Cleanup lcov
            COMMAND ${LCOV_PATH} ${LCOV_ARGS} --gcov-tool ${GCOV_PATH} -d ${PROJECT_BINARY_DIR} -b ${PROJECT_SOURCE_DIR} --zerocounters

            # Create output directory
            COMMAND ${CMAKE_COMMAND} -E make_directory ${LCOV_OUTPUT_DIR}

            # Create baseline to make sure untouched files show up in the report
            COMMAND ${LCOV_PATH} ${LCOV_ARGS} --gcov-tool ${GCOV_PATH} -d ${PROJECT_BINARY_DIR} -b ${PROJECT_SOURCE_DIR} -c -i -o ${LCOV_OUTPUT_DIR}/${LCOV_CASE_NAME}_base.info

            BYPRODUCTS ${LCOV_OUTPUT_DIR}
            VERBATIM
            COMMENT "Initializing lcov reports"
        )

        # Generate report
        add_custom_command(TARGET coverage POST_BUILD

            # Capturing lcov counters and generating report
            COMMAND ${LCOV_PATH} ${LCOV_ARGS} --gcov-tool ${GCOV_PATH} -d ${PROJECT_BINARY_DIR} -b ${PROJECT_SOURCE_DIR} -c -o ${LCOV_OUTPUT_DIR}/${LCOV_CASE_NAME}_test.info

            # Add baseline counters
            COMMAND ${LCOV_PATH} ${LCOV_ARGS} --gcov-tool ${GCOV_PATH} -a ${LCOV_OUTPUT_DIR}/${LCOV_CASE_NAME}_base.info -a ${LCOV_OUTPUT_DIR}/${LCOV_CASE_NAME}_test.info -o ${LCOV_OUTPUT_DIR}/${LCOV_CASE_NAME}_total.info

            # Filter collected data to final coverage report
            COMMAND ${LCOV_PATH} ${LCOV_ARGS} --gcov-tool ${GCOV_PATH} --remove ${LCOV_OUTPUT_DIR}/${LCOV_CASE_NAME}_total.info ${OVERALL_EXCLUDES} -o ${LCOV_OUTPUT_DIR}/${LCOV_CASE_NAME}.info

            # Generate HTML output
            COMMAND ${GENHTML_PATH} ${GENHTML_ARGS} -o ${LCOV_OUTPUT_DIR} ${LCOV_OUTPUT_DIR}/${LCOV_CASE_NAME}.info

            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
            VERBATIM
            COMMENT "Processing code coverage counters and generating report."
        )

        # Show path of the reports
        add_custom_command(TARGET coverage POST_BUILD
            COMMAND ;
            COMMENT "lcov coverage report saved in \"${LCOV_OUTPUT_DIR}/${LCOV_CASE_NAME}.info\"."
        )

        add_custom_command(TARGET coverage POST_BUILD
            COMMAND ;
            COMMENT "lcov HTML coverage report saved in \"${LCOV_OUTPUT_DIR}/index.html\"."
        )
    endif()
endfunction()

# ------------------------------------------------------------------------------------#
# Initialize library target
# ------------------------------------------------------------------------------------#
set(BITPIT_LIBRARY ${PROJECT_NAME} CACHE INTERNAL "bitpit library name" FORCE)

add_library(${BITPIT_LIBRARY})

if(MSVC)
    target_compile_features(${BITPIT_LIBRARY} PUBLIC cxx_std_11)
    set_target_properties(${BITPIT_LIBRARY} PROPERTIES CXX_STANDARD 11)
else()
    target_compile_features(${BITPIT_LIBRARY} PUBLIC cxx_std_14)
    set_target_properties(${BITPIT_LIBRARY} PROPERTIES CXX_STANDARD 14)
endif()

set_target_properties(${BITPIT_LIBRARY} PROPERTIES CXX_STANDARD_REQUIRED ON)

if(BITPIT_ENABLE_MPI)
    set_target_properties(${BITPIT_LIBRARY} PROPERTIES DEBUG_POSTFIX "_MPI_D")
    set_target_properties(${BITPIT_LIBRARY} PROPERTIES RELWITHDEBINFO_POSTFIX "_MPI_D")
    set_target_properties(${BITPIT_LIBRARY} PROPERTIES RELEASE_POSTFIX "_MPI")
else()
    set_target_properties(${BITPIT_LIBRARY} PROPERTIES DEBUG_POSTFIX "_D")
    set_target_properties(${BITPIT_LIBRARY} PROPERTIES RELWITHDEBINFO_POSTFIX "_D")
endif()

# ------------------------------------------------------------------------------------#
# Link Time Optimization (LTO)
# ------------------------------------------------------------------------------------#

# Detect if LTO can be enabled
if(DEFINED BITPIT_LTO_STRATEGY)
    if(NOT "${BITPIT_LTO_STRATEGY}" STREQUAL "Disabled")
        include(CheckIPOSupported)

        if("${BITPIT_LTO_STRATEGY}" STREQUAL "Enabled")
            if(NOT BUILD_SHARED_LIBS)
                message(FATAL_ERROR "LTO can be enabled only when building a shared library.")
            endif()

            if(MSYS OR MINGW)
                message(FATAL_ERROR "LTO in not currently supported on MinGW/MSYS.")
            endif()

            check_ipo_supported()
            set(ENABLE_LTO TRUE)
        elseif("${BITPIT_LTO_STRATEGY}" STREQUAL "Auto")
            if(NOT MSYS AND NOT MINGW AND NOT MSVC)
                if(BUILD_SHARED_LIBS)
                    check_ipo_supported(RESULT LTO_SUPPORTED)

                    if(${LTO_SUPPORTED} AND "${CMAKE_BUILD_TYPE}" STREQUAL "Release")
                        set(ENABLE_LTO TRUE)
                    else()
                        set(ENABLE_LTO FALSE)
                    endif()
                else()
                    set(ENABLE_LTO FALSE)
                endif()
            else()
                set(ENABLE_LTO FALSE)
            endif()
        else()
            set(ENABLE_LTO FALSE)
        endif()
    else()
        set(ENABLE_LTO FALSE)
    endif()
else()
    set(ENABLE_LTO FALSE)
endif()

# Set LTO property for the specified target
function(set_lto_property TARGET_NAME)
    if(${ENABLE_LTO})
        set_target_properties(${TARGET_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()
endfunction()

# Set LTO property for the main library
#
# Here LTO property is set only for the main library, each module may need
# to set its own LTO property.
set_lto_property(${BITPIT_LIBRARY})

# ------------------------------------------------------------------------------------#
# Sanitize
# ------------------------------------------------------------------------------------#
set(SANITIZE_MODULES "undefined;address")

# ------------------------------------------------------------------------------------#
# Customized build types
# ------------------------------------------------------------------------------------#

# Set build type variable
if(NOT CMAKE_BUILD_TYPE)
    if(MSVC)
        set(SUPPORTED_BUILD_TYPE "Coverage" "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
    else()
        set(SUPPORTED_BUILD_TYPE "Coverage" "Debug" "Release" "RelWithDebInfo" "RelWithSanitize" "MinSizeRel" "GNUProfiling" "ScalascaProfiling")
    endif()

    # Set default build type to Debug
    string(REPLACE ";" " " SUPPORTED_BUILD_TYPE_STRING "${SUPPORTED_BUILD_TYPE}")
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build, options are: ${SUPPORTED_BUILD_TYPE_STRING}." FORCE)
endif()

# ------------------------------------------------------------------------------------#
# Version
# ------------------------------------------------------------------------------------#

# Version regex expressions
set(VERSION_MATCHER_REGEX "^${CMAKE_PROJECT_NAME}-((([0-9]*)(\\.([0-9]*))?(\\.([0-9]*)))(-(.*))?)$")
set(VERSION_TAG_MATCHER_REGEX "^${CMAKE_PROJECT_NAME}-[^-]*(-(.*))?$")

# Initialize version string
set(BITPIT_VERSION_STRING "${CMAKE_PROJECT_NAME}-${BITPIT_VERSION}")

string(REGEX MATCH "${VERSION_MATCHER_REGEX}" VERSION_STRING_VALID "${BITPIT_VERSION_STRING}")

if("${VERSION_STRING_VALID}" STREQUAL "")
    message(FATAL_ERROR "Version string \"${BITPIT_VERSION_STRING}\" is not a valid version.")
endif()

# Initialize verbose version string
set(BITPIT_VERBOSE_VERSION_STRING "${BITPIT_VERSION_STRING}")

get_git_head_revision(BITPIT_GIT_REFERENCE BITPIT_GIT_HASH)

if(NOT "${BITPIT_GIT_REFERENCE}" STREQUAL "GITDIR-NOTFOUND")
    # Add git information
    #
    # Information about git are added only if git tag doens't match the version
    # string.
    git_get_exact_tag(BITPIT_GIT_TAG)

    if(NOT "${BITPIT_GIT_TAG}" STREQUAL "${BITPIT_VERBOSE_VERSION_STRING}")
        git_get_branch(BITPIT_GIT_BRANCH)

        if(NOT "${BITPIT_GIT_BRANCH}" MATCHES "NOTFOUND$")
            set(BITPIT_VERBOSE_VERSION_STRING "${BITPIT_VERBOSE_VERSION_STRING}-${BITPIT_GIT_BRANCH}")
        else()
            set(BITPIT_VERBOSE_VERSION_STRING "${BITPIT_VERBOSE_VERSION_STRING}-detached")
        endif()

        git_get_short_hash(BITPIT_GIT_SHORT_HASH)
        set(BITPIT_VERBOSE_VERSION_STRING "${BITPIT_VERBOSE_VERSION_STRING}-${BITPIT_GIT_SHORT_HASH}")
    endif()

    # Report if there are uncommited changes
    #
    # It's not possible to automatically re-run cmake (and thus update the
    # version) if there are uncommited changes. The label will be reliably
    # added only if, after changing the code, cmake is run manually.
    set(BITPIT_VERBOSE_VERSION_DIRTY_TAG "dirty" CACHE STRING "Tag that will be added to the verbose version string when there are uncommited changes")
    mark_as_advanced(BITPIT_VERBOSE_VERSION_DIRTY_TAG)

    git_local_changes(BITPIT_GIT_STATUS)

    if(BITPIT_GIT_STATUS STREQUAL "DIRTY")
        set(BITPIT_VERBOSE_VERSION_STRING "${BITPIT_VERBOSE_VERSION_STRING}-${BITPIT_VERBOSE_VERSION_DIRTY_LABEL}")
    endif()
endif()

string(REGEX MATCH "${VERSION_MATCHER_REGEX}" VERSION_STRING_VALID "${BITPIT_VERBOSE_VERSION_STRING}")

if("${VERSION_STRING_VALID}" STREQUAL "")
    message(FATAL_ERROR "Version string \"${BITPIT_VERBOSE_VERSION_STRING}\" is not a valid verbose version.")
endif()

# Extract version information
string(REGEX REPLACE "${VERSION_MATCHER_REGEX}" "\\1" BITPIT_VERBOSE_VERSION "${BITPIT_VERBOSE_VERSION_STRING}")

string(REGEX REPLACE "${VERSION_MATCHER_REGEX}" "\\1" BITPIT_VERSION "${BITPIT_VERSION_STRING}")
string(REGEX REPLACE "${VERSION_MATCHER_REGEX}" "\\3" BITPIT_MAJOR_VERSION "${BITPIT_VERSION_STRING}")
string(REGEX REPLACE "${VERSION_MATCHER_REGEX}" "\\5" BITPIT_MINOR_VERSION "${BITPIT_VERSION_STRING}")
string(REGEX REPLACE "${VERSION_MATCHER_REGEX}" "\\7" BITPIT_PATCH_VERSION "${BITPIT_VERSION_STRING}")

string(REGEX MATCH "${VERSION_TAG_MATCHER_REGEX}" dummy "${BITPIT_VERSION_STRING}")

if(NOT CMAKE_MATCH_2 STREQUAL "")
    string(REGEX REPLACE "${VERSION_TAG_MATCHER_REGEX}" "\\2" BITPIT_TAG_VERSION "${BITPIT_VERSION_STRING}")
else()
    set(BITPIT_TAG_VERSION "")
endif()

# Set library version
set_target_properties(${BITPIT_LIBRARY} PROPERTIES VERSION "${BITPIT_MAJOR_VERSION}.${BITPIT_MINOR_VERSION}.${BITPIT_PATCH_VERSION}"
    SOVERSION "${BITPIT_MAJOR_VERSION}")

# ------------------------------------------------------------------------------------#
# List of all avilable modules
# ------------------------------------------------------------------------------------#
set(BITPIT_MODULE_LIST "")
list(APPEND BITPIT_MODULE_LIST "common")
list(APPEND BITPIT_MODULE_LIST "operators")
list(APPEND BITPIT_MODULE_LIST "containers")
list(APPEND BITPIT_MODULE_LIST "IO")
list(APPEND BITPIT_MODULE_LIST "communications")
list(APPEND BITPIT_MODULE_LIST "LA")
list(APPEND BITPIT_MODULE_LIST "SA")
list(APPEND BITPIT_MODULE_LIST "CG")
list(APPEND BITPIT_MODULE_LIST "PABLO")
list(APPEND BITPIT_MODULE_LIST "patchkernel")
list(APPEND BITPIT_MODULE_LIST "pointcloud")
list(APPEND BITPIT_MODULE_LIST "lineunstructured")
list(APPEND BITPIT_MODULE_LIST "surfunstructured")
list(APPEND BITPIT_MODULE_LIST "volcartesian")
list(APPEND BITPIT_MODULE_LIST "voloctree")
list(APPEND BITPIT_MODULE_LIST "volunstructured")
list(APPEND BITPIT_MODULE_LIST "discretization")
list(APPEND BITPIT_MODULE_LIST "RBF")
list(APPEND BITPIT_MODULE_LIST "levelset")
list(APPEND BITPIT_MODULE_LIST "POD")

# ------------------------------------------------------------------------------------#
# Modules dependecies
# ------------------------------------------------------------------------------------#
set(OPERATORS_DEPS "common")
set(CONTAINERS_DEPS "common;operators")
set(IO_DEPS "common;operators;containers")
set(COMMUNICATIONS_DEPS "common;containers;IO")
set(LA_DEPS "common;containers;operators;IO")
set(SA_DEPS "common;operators")
set(CG_DEPS "common;operators;SA")
set(PABLO_DEPS "common;containers;IO;communications")
set(PATCHKERNEL_DEPS "common;operators;containers;IO;communications;CG")
set(POINTCLOUD_DEPS "common;patchkernel;IO")
set(LINEUNSTRUCTURED_DEPS "common;patchkernel;IO")
set(SURFUNSTRUCTURED_DEPS "common;patchkernel;lineunstructured;IO")
set(VOLCARTESIAN_DEPS "common;patchkernel")
set(VOLOCTREE_DEPS "common;PABLO;patchkernel")
set(VOLUNSTRUCTURED_DEPS "common;patchkernel")
set(RBF_DEPS "operators;IO")
set(DISCRETIZATION_DEPS "common;containers;LA;patchkernel")
set(LEVELSET_DEPS "common;communications;SA;CG;surfunstructured;voloctree;volcartesian;volunstructured;IO")
set(POD_DEPS "IO;common;containers;voloctree")

# ------------------------------------------------------------------------------------#
# Modules external dependecies
# ------------------------------------------------------------------------------------#
if(BITPIT_ENABLE_MPI)
    set(COMMON_EXTERNAL_DEPS "MPI")
endif()

set(OPERATORS_EXTERNAL_DEPS "")
set(CONTAINERS_EXTERNAL_DEPS "")
set(IO_EXTERNAL_DEPS "LibXml2;RapidJSON")
set(COMMUNICATIONS_EXTERNAL_DEPS "MPI")
set(LA_EXTERNAL_DEPS "PETSc")
set(SA_EXTERNAL_DEPS "")
set(CG_EXTERNAL_DEPS "")
set(PABLO_EXTERNAL_DEPS "")

if(BITPIT_ENABLE_MPI)
    set(PATCHKERNEL_EXTERNAL_DEPS "METIS")
else()
    set(PATCHKERNEL_EXTERNAL_DEPS "")
endif()

set(SURFUNSTRUCTURED_EXTERNAL_DEPS "")
set(VOLCARTESIAN_EXTERNAL_DEPS "")
set(VOLOCTREE_EXTERNAL_DEPS "")
set(VOLUNSTRUCTURED_EXTERNAL_DEPS "")
set(RBF_EXTERNAL_DEPS "LAPACKE")
set(DISCRETIZATION_EXTERNAL_DEPS "CBLAS;LAPACKE")
set(LEVELSET_EXTERNAL_DEPS "")
set(POD_EXTERNAL_DEPS "LAPACKE")

# ------------------------------------------------------------------------------------#
# Experimental/deprecated features
# ------------------------------------------------------------------------------------#

# ------------------------------------------------------------------------------------#
# Internal dependencies
# ------------------------------------------------------------------------------------#

# List of MPI-only modules
set(MPI_ONLY_MODULE_LIST "")

foreach(MODULE_NAME IN LISTS BITPIT_MODULE_LIST)
    string(TOUPPER ${MODULE_NAME} UPPER_MODULE_NAME)
    list(FIND ${UPPER_MODULE_NAME}_EXTERNAL_DEPS "MPI" _mpi_dep_index)

    if(${_mpi_dep_index} GREATER -1)
        LIST(APPEND MPI_ONLY_MODULE_LIST ${MODULE_NAME})
    endif()
endforeach()

# If building without MPI remove dependecies on MPI-only modules
if(NOT BITPIT_ENABLE_MPI)
    if(NOT ${MPI_ONLY_MODULE_LIST} STREQUAL "")
        foreach(MODULE_NAME IN LISTS BITPIT_MODULE_LIST)
            string(TOUPPER ${MODULE_NAME} UPPER_MODULE_NAME)
            list(LENGTH ${UPPER_MODULE_NAME}_DEPS DEPENDENCY_COUNT)

            if(${DEPENDENCY_COUNT} GREATER 0)
                list(REMOVE_ITEM ${UPPER_MODULE_NAME}_DEPS ${MPI_ONLY_MODULE_LIST})
            endif()

            unset(UPPER_MODULE_NAME)
        endforeach()
    endif()
endif()

# Initialize module state
foreach(MODULE_NAME IN LISTS BITPIT_MODULE_LIST)
    if((NOT BITPIT_ENABLE_MPI) AND(";${MPI_ONLY_MODULE_LIST};" MATCHES ";${MODULE_NAME};"))
        hideModule(${MODULE_NAME})
    else()
        enableModule(${MODULE_NAME} TRUE FALSE)
    endif()
endforeach()

set(CHECK_DEPS 1)

while(${CHECK_DEPS} EQUAL 1)
    set(CHECK_DEPS 0)

    foreach(MODULE_NAME IN LISTS BITPIT_MODULE_LIST)
        isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)

        if(MODULE_ENABLED)
            string(TOUPPER ${MODULE_NAME} UPPER_MODULE_NAME)

            foreach(DEPENDENCY_NAME IN LISTS ${UPPER_MODULE_NAME}_DEPS)
                isModuleEnabled(${DEPENDENCY_NAME} DEPENDENCY_ENABLED)

                if(NOT DEPENDENCY_ENABLED)
                    set(CHECK_DEPS 1)
                    enableModule(${DEPENDENCY_NAME} TRUE TRUE)
                    message("Module '${DEPENDENCY_NAME}' has been enabled because it's needed by module '${MODULE_NAME}'")
                endif()
            endforeach()
        endif()

        unset(MODULE_ENABLE_FLAG)
    endforeach()
endwhile()

unset(CHECK_DEPS)

# The modules 'common' and 'operators' has to be always built
getModuleEnableName("common" COMMON_MODULE_ENABLED)
unset(${COMMON_MODULE_ENABLED} CACHE)
set(${COMMON_MODULE_ENABLED} ON CACHE INTERNAL "Request building common module" FORCE)

getModuleEnableName("operators" OPERATORS_MODULE_ENABLED)
unset(${OPERATORS_MODULE_ENABLED} CACHE)
set(${OPERATORS_MODULE_ENABLED} ON CACHE INTERNAL "Request building operators module" FORCE)

# ------------------------------------------------------------------------------------#
# External dependecies
# ------------------------------------------------------------------------------------#

# List of dependencies
set(CBLAS_DEPS "BLAS")
set(LAPACK_DEPS "BLAS")
set(LAPACKE_DEPS "LAPACK")

if(BITPIT_ENABLE_MPI)
    set(PETSC_DEPS "MPI;LAPACK")
else()
    set(PETSC_DEPS "LAPACK")
endif()

set(EXTERNAL_DEPS_PROCESS_LIST "")

foreach(MODULE_NAME IN LISTS BITPIT_MODULE_LIST)
    isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)

    if(MODULE_ENABLED)
        string(TOUPPER ${MODULE_NAME} UPPER_MODULE_NAME)
        LIST(APPEND EXTERNAL_DEPS_PROCESS_LIST ${${UPPER_MODULE_NAME}_EXTERNAL_DEPS})
    endif()
endforeach()

set(EXTERNAL_DEPS "")
set(PROCESS_EXTERNAL_DEPS 1)

while(${PROCESS_EXTERNAL_DEPS} GREATER 0)
    set(EXTERNAL_DEPS_PROCESS_NEXT_LIST "")

    foreach(DEP_NAME IN LISTS EXTERNAL_DEPS_PROCESS_LIST)
        LIST(APPEND EXTERNAL_DEPS ${DEP_NAME})
        string(TOUPPER ${DEP_NAME} UPPER_DEP_NAME)

        if(DEFINED ${UPPER_DEP_NAME}_DEPS)
            LIST(APPEND EXTERNAL_DEPS_PROCESS_NEXT_LIST ${${UPPER_DEP_NAME}_DEPS})
        endif()
    endforeach()

    set(EXTERNAL_DEPS_PROCESS_LIST ${EXTERNAL_DEPS_PROCESS_NEXT_LIST})
    unset(EXTERNAL_DEPS_PROCESS_NEXT_LIST)
    list(LENGTH EXTERNAL_DEPS_PROCESS_LIST PROCESS_EXTERNAL_DEPS)
endwhile()

unset(PROCESS_EXTERNAL_DEPS)

list(REVERSE EXTERNAL_DEPS)
list(REMOVE_DUPLICATES EXTERNAL_DEPS)

# Enable additional languages
if(MSVC)
    enable_language("C")
else()
    set(BLAS_LANGUAGES "Fortran")
    set(CBLAS_LANGUAGES "C")
    set(LAPACK_LANGUAGES "Fortran")
    set(LAPACKE_LANGUAGES "C")
    set(METIS_LANGUAGES "C")
    set(PETSC_LANGUAGES "Fortran;C")

    foreach(DEP_NAME IN LISTS EXTERNAL_DEPS)
        if(DEFINED ${DEP_NAME}_LANGUAGES)
            string(TOUPPER ${DEP_NAME} UPPER_DEP_NAME)

            foreach(LANGUAGE IN LISTS ${UPPER_DEP_NAME}_LANGUAGES)
                enable_language(${LANGUAGE})
            endforeach()
        endif()
    endforeach()
endif()

# Add external dependencies
#
# Low-level dependencies are processed first.
set(BITPIT_EXTERNAL_DEPENDENCIES "")
set(BITPIT_EXTERNAL_VARIABLES_LIBRARIES "")
set(BITPIT_EXTERNAL_VARIABLES_INCLUDE_DIRS "")

list(FIND EXTERNAL_DEPS "MPI" _MPI_index)

if(${_MPI_index} GREATER -1)
    find_package(MPI)

    target_compile_definitions(${BITPIT_LIBRARY} PUBLIC "BITPIT_ENABLE_MPI=1")

    if(MPI_CXX_COMPILE_FLAGS)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
    endif()

    if(MPI_CXX_LINK_FLAGS)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")
    endif()

    list(INSERT BITPIT_EXTERNAL_DEPENDENCIES 0 "MPI")
    list(INSERT BITPIT_EXTERNAL_VARIABLES_LIBRARIES 0 "MPI_CXX_LIBRARIES")

    if("Fortran" IN_LIST BITPIT_LANGUAGES)
        list(INSERT BITPIT_EXTERNAL_VARIABLES_LIBRARIES 0 "MPI_Fortran_LIBRARIES")
    endif()

    list(INSERT BITPIT_EXTERNAL_VARIABLES_INCLUDE_DIRS 0 "MPI_CXX_INCLUDE_PATH")
else()
    target_compile_definitions(${BITPIT_LIBRARY} PUBLIC "BITPIT_ENABLE_MPI=0")
endif()

unset(_MPI_index)

if(MSVC)
    # download directly libxml2, rapidjson, mkl blas, mkl lapack and mkl lapacke
    # through nuget system and attach it to current VS project
    include("${CMAKE_CURRENT_LIST_DIR}/external/CMakeNuGetTools/NuGetTools.cmake")

    # find the nuget.exe executable
    find_program(NUGET_COMMAND nuget)

    if(NOT NUGET_COMMAND)
        message(FATAL "Cannot find nuget command line tool that should come with VS!")
    endif()

    # Call this once before any other nuget_* calls.
    nuget_initialize()

    # NuGet will install .
    nuget_add_dependencies(
        PACKAGE libxml2 VERSION 2.7.8.7 CMAKE_PREFIX_PATHS build/native;build/native/lib/v110/x64/Release/dynamic/cdecl
        PACKAGE libiconv VERSION 1.14.0.11 CMAKE_PREFIX_PATHS build/native
        PACKAGE zlib_native.redist VERSION 1.2.11 CMAKE_PREFIX_PATHS build/native
        PACKAGE tencent.rapidjson VERSION 1.1.1 CMAKE_PREFIX_PATHS lib/native
        PACKAGE intelmkl.devel.win-x64 VERSION 2022.0.0.115 CMAKE_PREFIX_PATHS lib/native;lib/native/win-x64
    )

    # # add C library missing GNU libiconv header
    include_directories("${NUGET_PACKAGES_DIR}/libiconv/build/native/include")

    # #TEMPORARY
    # # set temporary exception to include intelmkl include dirs
    set(INTELMKL_INCLUDE_DIRS "${NUGET_PACKAGES_DIR}/intelmkl.devel.win-x64/lib/native/include")

    # create ad-hoc CMakeConfig for tencent.rapidjson (they do not provide it in nuget package)
    set(CONFIGCMAKEFILENAME "${NUGET_PACKAGES_DIR}/tencent.rapidjson/lib/native/RapidJSONConfig.cmake")
    file(WRITE ${CONFIGCMAKEFILENAME} "get_filename_component(RAPIDJSON_CMAKE_DIR \"\${CMAKE_CURRENT_LIST_FILE}\" PATH)\n")
    file(APPEND ${CONFIGCMAKEFILENAME} "set(RAPIDJSON_INCLUDE_DIRS \"\${RAPIDJSON_CMAKE_DIR}/include\")\n")
    file(APPEND ${CONFIGCMAKEFILENAME} "message(STATUS \"RapidJSON found. Headers: \${RAPIDJSON_INCLUDE_DIRS}\")")

    set(CONFIGCMAKEFILEVERSIONNAME "${NUGET_PACKAGES_DIR}/tencent.rapidjson/lib/native/RapidJSONConfigVersion.cmake")
    file(WRITE ${CONFIGCMAKEFILEVERSIONNAME} "SET(PACKAGE_VERSION \"1.1.0\")\n\n")
    file(APPEND ${CONFIGCMAKEFILEVERSIONNAME} "IF (PACKAGE_FIND_VERSION VERSION_EQUAL PACKAGE_VERSION)\n")
    file(APPEND ${CONFIGCMAKEFILEVERSIONNAME} "    SET(PACKAGE_VERSION_EXACT \"true\")\n")
    file(APPEND ${CONFIGCMAKEFILEVERSIONNAME} "ENDIF (PACKAGE_FIND_VERSION VERSION_EQUAL PACKAGE_VERSION)\n")
    file(APPEND ${CONFIGCMAKEFILEVERSIONNAME} "IF (NOT PACKAGE_FIND_VERSION VERSION_GREATER PACKAGE_VERSION)\n")
    file(APPEND ${CONFIGCMAKEFILEVERSIONNAME} "    SET(PACKAGE_VERSION_COMPATIBLE \"true\")\n")
    file(APPEND ${CONFIGCMAKEFILEVERSIONNAME} "ELSE (NOT PACKAGE_FIND_VERSION VERSION_GREATER PACKAGE_VERSION)\n")
    file(APPEND ${CONFIGCMAKEFILEVERSIONNAME} "    SET(PACKAGE_VERSION_UNSUITABLE \"true\")\n")
    file(APPEND ${CONFIGCMAKEFILEVERSIONNAME} "ENDIF (NOT PACKAGE_FIND_VERSION VERSION_GREATER PACKAGE_VERSION)")

    list(APPEND CMAKE_MODULE_PATH "${NUGET_PACKAGES_DIR}/tencent.rapidjson/lib/native/")
endif()

list(FIND EXTERNAL_DEPS "BLAS" _BLAS_index)

if(${_BLAS_index} GREATER -1)
    set(BLAS_VENDOR "All" CACHE STRING "If set, checks only the specified vendor. If not set, checks all the possibilities")
    mark_as_advanced(BLAS_VENDOR)
    set(BLA_VENDOR "${BLAS_VENDOR}")

    find_package(BLAS REQUIRED)

    list(INSERT BITPIT_EXTERNAL_DEPENDENCIES 0 "BLAS")
    list(INSERT BITPIT_EXTERNAL_VARIABLES_LIBRARIES 0 "BLAS_LIBRARIES")

    if(MSVC)
        list(FIND BITPIT_EXTERNAL_VARIABLES_INCLUDE_DIRS "INTELMKL_INCLUDE_DIRS" _Intelmkl_index)

        if(${_Intelmkl_index} EQUAL -1)
            list(INSERT BITPIT_EXTERNAL_VARIABLES_INCLUDE_DIRS 0 "INTELMKL_INCLUDE_DIRS")
        endif()

        unset(_Intelmkl_index)
    else()
        list(INSERT BITPIT_EXTERNAL_VARIABLES_INCLUDE_DIRS 0 "BLAS_INCLUDE_DIRS")
    endif()
endif()

unset(_BLAS_index)

list(FIND EXTERNAL_DEPS "CBLAS" _CBLAS_index)

if(${_CBLAS_index} GREATER -1)
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/external/CBLAS/cmake")

    find_package(CBLAS REQUIRED)

    list(INSERT BITPIT_EXTERNAL_DEPENDENCIES 0 "CBLAS")
    list(INSERT BITPIT_EXTERNAL_VARIABLES_LIBRARIES 0 "CBLAS_LIBRARIES")

    if(MSVC)
        list(FIND BITPIT_EXTERNAL_VARIABLES_INCLUDE_DIRS "INTELMKL_INCLUDE_DIRS" _Intelmkl_index)

        if(${_Intelmkl_index} EQUAL -1)
            list(INSERT BITPIT_EXTERNAL_VARIABLES_INCLUDE_DIRS 0 "INTELMKL_INCLUDE_DIRS")
        endif()

        unset(_Intelmkl_index)
    else()
        list(INSERT BITPIT_EXTERNAL_VARIABLES_INCLUDE_DIRS 0 "CBLAS_INCLUDE_DIRS")
    endif()
endif()

unset(_CBLAS_index)

list(FIND EXTERNAL_DEPS "LAPACK" _LAPACK_index)

if(${_LAPACK_index} GREATER -1)
    find_package(LAPACK REQUIRED)

    list(INSERT BITPIT_EXTERNAL_DEPENDENCIES 0 "LAPACK")
    list(INSERT BITPIT_EXTERNAL_VARIABLES_LIBRARIES 0 "LAPACK_LIBRARIES")

    if(MSVC)
        list(FIND BITPIT_EXTERNAL_VARIABLES_INCLUDE_DIRS "INTELMKL_INCLUDE_DIRS" _Intelmkl_index)

        if(${_Intelmkl_index} EQUAL -1)
            list(INSERT BITPIT_EXTERNAL_VARIABLES_INCLUDE_DIRS 0 "INTELMKL_INCLUDE_DIRS")
        endif()

        unset(_Intelmkl_index)
    else()
        list(INSERT BITPIT_EXTERNAL_VARIABLES_INCLUDE_DIRS 0 "LAPACK_INCLUDE_DIRS")
    endif()
endif()

unset(_LAPACK_index)

list(FIND EXTERNAL_DEPS "LAPACKE" _LAPACKE_index)

if(${_LAPACKE_index} GREATER -1)
    set(LAPACKE_STATIC 0 CACHE BOOL "Require LAPACKE static libraries")
    mark_as_advanced(LAPACKE_STATIC)

    list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/external/LAPACKE/cmake")

    find_package(LAPACKE REQUIRED)

    set(LAPACKE_PRIVATE_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/external/LAPACKE/include")

    list(INSERT BITPIT_EXTERNAL_DEPENDENCIES 0 "LAPACKE")
    list(INSERT BITPIT_EXTERNAL_VARIABLES_LIBRARIES 0 "LAPACKE_LIBRARIES")
    list(INSERT BITPIT_EXTERNAL_VARIABLES_INCLUDE_DIRS 0 "LAPACKE_PRIVATE_INCLUDE_DIR")

    if(MSVC)
        list(FIND BITPIT_EXTERNAL_VARIABLES_INCLUDE_DIRS "INTELMKL_INCLUDE_DIRS" _Intelmkl_index)

        if(${_Intelmkl_index} EQUAL -1)
            list(INSERT BITPIT_EXTERNAL_VARIABLES_INCLUDE_DIRS 0 "INTELMKL_INCLUDE_DIRS")
        endif()

        unset(_Intelmkl_index)
    else()
        list(INSERT BITPIT_EXTERNAL_VARIABLES_INCLUDE_DIRS 0 "LAPACKE_INCLUDE_DIRS")
    endif()
endif()

unset(_LAPACKE_index)

list(FIND EXTERNAL_DEPS "LibXml2" _LibXml2_index)

if(${_LibXml2_index} GREATER -1)
    find_package(LibXml2 REQUIRED)

    list(INSERT BITPIT_EXTERNAL_DEPENDENCIES 0 "LibXml2")
    list(INSERT BITPIT_EXTERNAL_VARIABLES_LIBRARIES 0 "LIBXML2_LIBRARIES")
    list(INSERT BITPIT_EXTERNAL_VARIABLES_INCLUDE_DIRS 0 "LIBXML2_INCLUDE_DIR")
endif()

unset(_LibXml2_index)

list(FIND EXTERNAL_DEPS "RapidJSON" _RapidJSON_index)

if(${_RapidJSON_index} GREATER -1)
    find_package(RapidJSON 1.1.0 QUIET)
    mark_as_advanced(RapidJSON_DIR)

    if(RapidJSON_FOUND)
        # RapidJSON library has been found
        target_compile_definitions(${BITPIT_LIBRARY} PUBLIC "BITPIT_ENABLE_RAPIDJSON=1")

        list(INSERT BITPIT_EXTERNAL_DEPENDENCIES 0 "RapidJSON")
        list(INSERT BITPIT_EXTERNAL_VARIABLES_INCLUDE_DIRS 0 "RAPIDJSON_INCLUDE_DIRS")

        # Activate std::string RapidJSON compliancy
        target_compile_definitions(${BITPIT_LIBRARY} PUBLIC "RAPIDJSON_HAS_STDSTRING=1")
    else()
        # RapidJSON library has not been found
        target_compile_definitions(${BITPIT_LIBRARY} PUBLIC "BITPIT_ENABLE_RAPIDJSON=0")

        message(STATUS "RapidJSON library not found, JSON support will be disabled.")
    endif()
endif()

unset(_RapidJSON_index)

list(FIND EXTERNAL_DEPS "PETSc" _PETSc_index)

if(${_PETSc_index} GREATER -1)
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/external/PETSc/cmake")

    if(NOT PETSC_DIR AND DEFINED ENV{PETSC_DIR})
        set(DEFAULT_PETSC_DIR "$ENV{PETSC_DIR}")
    else()
        set(DEFAULT_PETSC_DIR "")
    endif()

    set(PETSC_DIR "${DEFAULT_PETSC_DIR}" CACHE PATH "Installation directory of PETSC library")

    if(NOT PETSC_ARCH AND DEFINED ENV{PETSC_ARCH})
        set(DEFAULT_PETSC_ARCH "$ENV{PETSC_ARCH}")
    else()
        set(DEFAULT_PETSC_ARCH "")
    endif()

    set(PETSC_ARCH "${DEFAULT_PETSC_ARCH}" CACHE STRING "Build architecture")

    find_package(PETSc REQUIRED)

    list(INSERT BITPIT_EXTERNAL_DEPENDENCIES 0 "PETSc")
    list(INSERT BITPIT_EXTERNAL_VARIABLES_LIBRARIES 0 "PETSC_LIBRARIES")
    list(INSERT BITPIT_EXTERNAL_VARIABLES_INCLUDE_DIRS 0 "PETSC_INCLUDES")
endif()

unset(_PETSc_index)

# list(FIND EXTERNAL_DEPS "BoostTest" _BoostTest_index)

# if(${_BoostTest_index} GREATER -1)
# set(BOOST_COMPONENTS "")

# if(${_BoostTest_index} GREATER -1)
# LIST(APPEND BOOST_COMPONENTS "unit_test_framework")
# endif()

# string(REPLACE ";" " " BOOST_COMPONENTS_STRING "${BOOST_COMPONENTS}")

# FIND_PACKAGE(Boost 1.40 COMPONENTS "${BOOST_COMPONENTS_STRING}" REQUIRED)

# list(INSERT BITPIT_EXTERNAL_DEPENDENCIES 0 "Boost")
# list(INSERT BITPIT_EXTERNAL_VARIABLES_LIBRARIES 0 "Boost_LIBRARIES")
# list(INSERT BITPIT_EXTERNAL_VARIABLES_INCLUDE_DIRS 0 "Boost_INCLUDE_DIR")
# endif()

# unset(_BoostTest_index)
list(FIND EXTERNAL_DEPS "METIS" _METIS_index)

if(${_METIS_index} GREATER -1)
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/external/METIS/cmake")

    find_package(METIS)

    if(METIS_FOUND)
        target_compile_definitions(${BITPIT_LIBRARY} PUBLIC "BITPIT_ENABLE_METIS=1")

        list(APPEND BITPIT_EXTERNAL_DEPENDENCIES "METIS")
        list(APPEND BITPIT_EXTERNAL_VARIABLES_LIBRARIES "METIS_LIBRARIES")
        list(APPEND BITPIT_EXTERNAL_VARIABLES_INCLUDE_DIRS "METIS_INCLUDE_DIRS")
    else()
        target_compile_definitions(${BITPIT_LIBRARY} PUBLIC "BITPIT_ENABLE_METIS=0")

        message(STATUS "METIS library not found, default patch partitioning support will be disabled.")
    endif()
endif()

unset(_METIS_index)

set(BITPIT_EXTERNAL_LIBRARIES "")

foreach(VARIABLE_NAME IN LISTS BITPIT_EXTERNAL_VARIABLES_LIBRARIES)
    list(APPEND BITPIT_EXTERNAL_LIBRARIES "${${VARIABLE_NAME}}")
endforeach()

set(BITPIT_EXTERNAL_INCLUDE_DIRS "")

foreach(VARIABLE_NAME IN LISTS BITPIT_EXTERNAL_VARIABLES_INCLUDE_DIRS)
    list(APPEND BITPIT_EXTERNAL_INCLUDE_DIRS "${${VARIABLE_NAME}}")
endforeach()

# Windows shared libraries have to be linked against their dependencies.
if((MSVC OR MINGW) AND BUILD_SHARED_LIBS)
    target_link_libraries(${BITPIT_LIBRARY} ${BITPIT_EXTERNAL_LIBRARIES})
endif()

# ------------------------------------------------------------------------------------#
# Compiler settings
# ------------------------------------------------------------------------------------#
string(TOLOWER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_LOWER)

# GNU profiling build type
SET(CMAKE_CXX_FLAGS_GNUPROFILING "-pg" CACHE STRING
    "Flags used by the C++ compiler during GNU profiling builds." FORCE)
MARK_AS_ADVANCED(CMAKE_CXX_FLAGS_GNUPROFILING)

SET(CMAKE_C_FLAGS_GNUPROFILING "-pg" CACHE STRING
    "Flags used by the C compiler during GNU profiling builds." FORCE)
MARK_AS_ADVANCED(CMAKE_C_FLAGS_GNUPROFILING)

SET(CMAKE_EXE_LINKER_FLAGS_GNUPROFILING "-pg" CACHE STRING
    "Flags used for linking binaries during GNU profiling builds." FORCE)
MARK_AS_ADVANCED(CMAKE_EXE_LINKER_FLAGS_GNUPROFILING)

SET(CMAKE_SHARED_LINKER_FLAGS_GNUPROFILING "-pg" CACHE STRING
    "Flags used by the shared libraries linker during GNU profiling builds." FORCE)
MARK_AS_ADVANCED(CMAKE_SHARED_LINKER_FLAGS_GNUPROFILING)

# Scalasca profiling build type
SET(CMAKE_CXX_FLAGS_SCALASCAPROFILING "-O2" CACHE STRING
    "Flags used by the C++ compiler during Scalasca profiling builds." FORCE)
MARK_AS_ADVANCED(CMAKE_CXX_FLAGS_SCALASCAPROFILING)

SET(CMAKE_C_FLAGS_SCALASCAPROFILING "-O2" CACHE STRING
    "Flags used by the C compiler during Scalasca builds." FORCE)
MARK_AS_ADVANCED(CMAKE_C_FLAGS_SCALASCAPROFILING)

SET(CMAKE_EXE_LINKER_FLAGS_SCALASCAPROFILING "" CACHE STRING
    "Flags used for linking binaries during Scalasca builds." FORCE)
MARK_AS_ADVANCED(CMAKE_EXE_LINKER_FLAGS_SCALASCAPROFILING)

SET(CMAKE_SHARED_LINKER_FLAGS_SCALASCAPROFILING "" CACHE STRING
    "Flags used by the shared libraries linker during Scalasca builds." FORCE)
MARK_AS_ADVANCED(CMAKE_SHARED_LINKER_FLAGS_SCALASCAPROFILING)

# Standard build types
if("C" IN_LIST BITPIT_LANGUAGES)
    set(C_SANITIZE_MODULE_FLAGS "")

    foreach(SANITIZE_MODULE IN LISTS SANITIZE_MODULES)
        check_c_compiler_flag(-fsanitize=${SANITIZE_MODULE} C_COMPILER_SUPPORTS_SANITIZE_MODULE)

        if(C_COMPILER_SUPPORTS_SANITIZE_MODULE)
            set(C_SANITIZE_MODULE_FLAGS "${C_SANITIZE_MODULE_FLAGS} -fsanitize=${SANITIZE_MODULE}")
        endif()
    endforeach()

    if(CMAKE_BUILD_TYPE_LOWER MATCHES "sanitize")
        if("${C_SANITIZE_MODULE_FLAGS}" STREQUAL "")
            message(FATAL_ERROR "The compiler doesn't support sanitize")
        endif()
    endif()
endif()

set(CXX_SANITIZE_MODULE_FLAGS "")

foreach(SANITIZE_MODULE IN LISTS SANITIZE_MODULES)
    check_cxx_compiler_flag(-fsanitize=${SANITIZE_MODULE} CXX_COMPILER_SUPPORTS_SANITIZE_MODULE)

    if(CXX_COMPILER_SUPPORTS_SANITIZE_MODULE)
        set(CXX_SANITIZE_MODULE_FLAGS "${CXX_SANITIZE_MODULE_FLAGS} -fsanitize=${SANITIZE_MODULE}")
    endif()
endforeach()

if(CMAKE_BUILD_TYPE_LOWER MATCHES "sanitize")
    if("${CXX_SANITIZE_MODULE_FLAGS}" STREQUAL "")
        message(FATAL_ERROR "The compiler doesn't support sanitize")
    endif()
endif()

# Standard build types
if(MSVC)
    set(CMAKE_C_FLAGS_RELWITHDEBINFO "/O2 /Zi /Zo")
    set(CMAKE_C_FLAGS_DEBUG "/Od /Zi")
    set(CMAKE_C_FLAGS_RELEASE "/O2")

    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "/O2 /Zi /Zo")
    set(CMAKE_CXX_FLAGS_DEBUG "/Od /Zi")
    set(CMAKE_CXX_FLAGS_RELEASE "/O2")

    # keep all warning always on - this is the policy of bitpit.
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W3")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3")
else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0")
    set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g")
    set(CMAKE_C_FLAGS_RELWITHSANITIZE "-O2 -g -fno-omit-frame-pointer -fno-common ${C_SANITIZE_MODULE_FLAGS}")
    set(CMAKE_C_FLAGS_DEBUG "-O0 -g")
    set(CMAKE_C_FLAGS_RELEASE "-O2")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
    set(CMAKE_CXX_FLAGS_RELWITHSANITIZE "-O2 -g -fno-omit-frame-pointer -fno-common ${CXX_SANITIZE_MODULE_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2")

    CHECK_CXX_COMPILER_FLAG("-Wall" CXX_SUPPORTS_WALL)

    if(CXX_SUPPORTS_WALL)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    endif()

    CHECK_CXX_COMPILER_FLAG("-Wextra" CXX_SUPPORTS_WEXTRA)

    if(CXX_SUPPORTS_WEXTRA)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
    endif()

    CHECK_CXX_COMPILER_FLAG("-Wsuggest-override" CXX_SUPPORTS_WSUGGEST_OVERRIDE)

    if(CXX_SUPPORTS_WSUGGEST_OVERRIDE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsuggest-override")
    endif()

    CHECK_C_COMPILER_FLAG("-Wall" C_SUPPORTS_WALL)

    if(C_SUPPORTS_WALL)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
    endif()

    CHECK_C_COMPILER_FLAG("-Wextra" C_SUPPORTS_WEXTRA)

    if(C_SUPPORTS_WEXTRA)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wextra")
    endif()

    set(CMAKE_EXE_LINKER_FLAGS_RELWITHSANITIZE "-fno-omit-frame-pointer -fno-common ${CXX_SANITIZE_MODULE_FLAGS}")
endif()

# Coverage build type
SET(CMAKE_C_FLAGS_COVERAGE "-g -O0 -fprofile-arcs -ftest-coverage" CACHE STRING
    "Flags used by the C compiler during coverage builds." FORCE)
MARK_AS_ADVANCED(CMAKE_C_FLAGS_COVERAGE)

SET(CMAKE_CXX_FLAGS_COVERAGE "-g -O0 -fprofile-arcs -ftest-coverage" CACHE STRING
    "Flags used by the C++ compiler during coverage builds." FORCE)
MARK_AS_ADVANCED(CMAKE_CXX_FLAGS_COVERAGE)

SET(CMAKE_EXE_LINKER_FLAGS_COVERAGE "-fprofile-arcs -ftest-coverage" CACHE STRING
    "Flags used for linking binaries during coverage builds." FORCE)
MARK_AS_ADVANCED(CMAKE_EXE_LINKER_FLAGS_COVERAGE)

SET(CMAKE_SHARED_LINKER_FLAGS_COVERAGE "-fprofile-arcs -ftest-coverage" CACHE STRING
    "Flags used by the shared libraries linker during coverage builds." FORCE)
MARK_AS_ADVANCED(CMAKE_SHARED_LINKER_FLAGS_COVERAGE)

# Set linker flags for Intel compilers
#
# If we are using an external library that uses Fortran, we need to add the
# Fortran core library to the linker flags. This is somewhat a workaround,
# because, for building C++ applications that call Fortran functions, Intel
# recommends to link the C++ program with the Fortran compiler (ifort) and
# passing it the flags "-cxxlib" and "-nofor_main". This would be difficult
# to achieve with CMake, so we link with the C++ linker and we pass it the
# Fortran core library.
if(CMAKE_C_COMPILER_ID STREQUAL "Intel" OR CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    if("Fortran" IN_LIST BITPIT_LANGUAGES)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lifcore")
    endif()
endif()

# Define a preprocessor macro to recognize debug builds
IF(CMAKE_BUILD_TYPE_LOWER MATCHES "debug" OR CMAKE_BUILD_TYPE_LOWER MATCHES "debinfo")
    target_compile_definitions(${BITPIT_LIBRARY} PUBLIC "BITPIT_ENABLE_DEBUG=1")
else()
    target_compile_definitions(${BITPIT_LIBRARY} PUBLIC "BITPIT_ENABLE_DEBUG=0")
endif()

IF(NOT CMAKE_BUILD_TYPE_LOWER MATCHES "debug")
    target_compile_definitions(${BITPIT_LIBRARY} PRIVATE "NDEBUG")
endif()

# Define an alias for building with scalasca
if(BITPIT_ENABLE_MPI)
    SET(C_FLAGS_INSTRUMENT "-instrument mpicxx")
    SET(CXX_FLAGS_INSTRUMENT "-instrument mpic")
else()
    SET(C_FLAGS_INSTRUMENT "")
    SET(CXX_FLAGS_INSTRUMENT "")
endif()

if(CMAKE_BUILD_TYPE_LOWER MATCHES "scalasca")
    file(WRITE scalasca_c_compiler
        "#!/bin/bash
scalasca ${C_FLAGS_INSTRUMENT} \"$@\"
"
    )

    file(WRITE scalasca_cxx_compiler
        "#!/bin/bash
scalasca ${C_FLAGS_INSTRUMENT} \"$@\"
"
    )

    file(INSTALL scalasca_cxx_compiler DESTINATION ${PROJECT_BINARY_DIR} PERMISSIONS OWNER_READ OWNER_EXECUTE)
    file(INSTALL scalasca_c_compiler DESTINATION ${PROJECT_BINARY_DIR} PERMISSIONS OWNER_READ OWNER_EXECUTE)
    file(REMOVE scalasca_cxx_compiler)
    file(REMOVE scalasca_c_compiler)

    set(CMAKE_CXX_COMPILER "${PROJECT_BINARY_DIR}/scalasca_cxx_compiler")
    set(CMAKE_C_COMPILER "${PROJECT_BINARY_DIR}/scalasca_c_compiler")
endif()

# Check the features supported by the compiler
include(CheckCXXSourceCompiles)

CHECK_CXX_SOURCE_COMPILES("int main() {__builtin_unreachable();}" BITPIT_ENABLE_BUILTIN_UNREACHABLE)

if(BITPIT_ENABLE_BUILTIN_UNREACHABLE)
    target_compile_definitions(${BITPIT_LIBRARY} PUBLIC "BITPIT_ENABLE_BUILTIN_UNREACHABLE")
endif()

# ------------------------------------------------------------------------------------#
# Source code formatting
# ------------------------------------------------------------------------------------#
SET(CLANG_FORMAT_TARGET "clang-format")
ClangFormatAddTarget(${CLANG_FORMAT_TARGET})

# ------------------------------------------------------------------------------------#
# Code coverage analysis
# ------------------------------------------------------------------------------------#
if(CMAKE_BUILD_TYPE STREQUAL "Coverage")
    if(NOT "${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
        message(FATAL_ERROR "Code coverage analysis requires gcc! Aborting...")
    endif()

    set(BITPIT_COVERAGE_GCOVR_REPORTS ON CACHE BOOL "If set, coverage analysis will generate gcov reports")
    mark_as_advanced(${BITPIT_COVERAGE_GCOVR_REPORTS})

    set(BITPIT_COVERAGE_LCOV_REPORTS ON CACHE BOOL "If set, coverage analysis will generate lcov reports")
    mark_as_advanced(${BITPIT_COVERAGE_LCOV_REPORTS})

    addCoverageAnalysis(${BITPIT_COVERAGE_GCOVR_REPORTS} ${BITPIT_COVERAGE_LCOV_REPORTS} "${PROJECT_SOURCE_DIR}/examples/*;${PROJECT_SOURCE_DIR}/test/*")
endif()

# ------------------------------------------------------------------------------------#
# Set programming languages
# ------------------------------------------------------------------------------------#
get_property(BITPIT_LANGUAGES GLOBAL PROPERTY ENABLED_LANGUAGES)

# ------------------------------------------------------------------------------------#
# Subdirectories
# ------------------------------------------------------------------------------------#

# Main library
add_subdirectory(src)

# CMake
add_subdirectory(cmake)

# Examples
add_subdirectory(examples)

# External
add_subdirectory(external)

# Docs
if(BITPIT_ENABLE_DOC)
    add_subdirectory(doc)
endif()

# Tests
if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(test)
endif()

# ------------------------------------------------------------------------------------#
# Installation
# ------------------------------------------------------------------------------------#
install(TARGETS ${BITPIT_LIBRARY}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_PROJECT_NAME}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_PROJECT_NAME}
)
