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

#------------------------------------------------------------------------------------#
# CMake settings
#------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 2.8)

project("bitpit")

#------------------------------------------------------------------------------------#
# Variables visible to the user
#------------------------------------------------------------------------------------#
set(ENABLE_MPI 0 CACHE BOOL "If set, the program is compiled with MPI support")
set(VERBOSE_MAKE 0 CACHE BOOL "Set appropriate compiler and cmake flags to enable verbose output from compilation")

#------------------------------------------------------------------------------------#
# Internal variables
#------------------------------------------------------------------------------------#
set(BUILD_SHARED_LIBS 0 CACHE INTERNAL "Build Shared Libraries")

#------------------------------------------------------------------------------------#
# Functions
#------------------------------------------------------------------------------------#
function(getModuleEnableName MODULE_NAME FLAG_NAME)
	string(TOUPPER ${MODULE_NAME} UPPER_MODULE_NAME)
	set(${FLAG_NAME} "BITPIT_MODULE_${UPPER_MODULE_NAME}" PARENT_SCOPE)
endfunction()

function(isModuleEnabled MODULE_NAME ENABLED)
	getModuleEnableName(${MODULE_NAME} ENABLED_VARIABLE)
	if (DEFINED ${ENABLED_VARIABLE})
		set(${ENABLED} ${${ENABLED_VARIABLE}} PARENT_SCOPE)
	else ()
		set(${ENABLED} "FALSE" PARENT_SCOPE)
	endif ()
endfunction()

function(areModulesEnabled MODULE_LIST ENABLED)
	foreach (MODULE_NAME IN LISTS MODULE_LIST)
		isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)
		if (NOT MODULE_ENABLED)
			set(${ENABLED} "FALSE" PARENT_SCOPE)
			return()
		endif()
	endforeach ()

	set(${ENABLED} "TRUE" PARENT_SCOPE)
endfunction()

function(enableModule MODULE_NAME ENABLED)
	set (EXTRA_ARGUMENTS ${ARGN})
	list(LENGTH EXTRA_ARGUMENTS EXTRA_ARGUMENT_COUNT)
	if (${EXTRA_ARGUMENT_COUNT} GREATER 0)
		list(GET EXTRA_ARGUMENTS 0 FORCED)
		if (FORCED)
			set(FORCE_FLAG "FORCE")
		endif()
	endif ()

	if (NOT DEFINED FORCE_FLAG)
		set(FORCE_FLAG "")
	endif ()

	getModuleEnableName(${MODULE_NAME} MODULE_ENABLE_FLAG)
	set(${MODULE_ENABLE_FLAG} ${ENABLED} CACHE BOOL "Request building ${MODULE_NAME} module" ${FORCE_FLAG})
	mark_as_advanced(${MODULE_ENABLE_FLAG})
endfunction()

function(hideModule MODULE_NAME)
	getModuleEnableName(${MODULE_NAME} MODULE_ENABLE_FLAG)
	if (DEFINED ${MODULE_ENABLE_FLAG})
		unset(${MODULE_ENABLE_FLAG} CACHE)
	endif ()
endfunction()

function(getTestName TEST_ENTRY TEST_NAME)
	string(REGEX REPLACE ":[^:]*$" "" TEST_ENTRY_STRIPPED ${TEST_ENTRY})
	set(${TEST_NAME} "${TEST_ENTRY_STRIPPED}" PARENT_SCOPE)
endfunction()

function(isTestParallel TEST_ENTRY PARALLEL_FLAG)
	if (TEST_ENTRY MATCHES "parallel")
		set(${PARALLEL_FLAG} "1" PARENT_SCOPE)
	else ()
		set(${PARALLEL_FLAG} "0" PARENT_SCOPE)
	endif ()
endfunction()

function(getTestNumProcs TEST_ENTRY DEFAULT_N_PROCS N_PROCS)
	isTestParallel("${TEST_ENTRY}" IS_TEST_PARALLEL)
	if (NOT IS_TEST_PARALLEL)
		set(${N_PROCS} "1" PARENT_SCOPE)
		return()
	endif ()

	string(REGEX MATCH "[^:]*$" TEST_ENTRY_MATCH ${TEST_ENTRY})
	if ("${TEST_ENTRY_MATCH}" STREQUAL "${TEST_ENTRY}")
		set(${N_PROCS} "${DEFAULT_N_PROCS}" PARENT_SCOPE)
		return()
	endif ()

	set(${N_PROCS} "${TEST_ENTRY_MATCH}" PARENT_SCOPE)
endfunction()

#------------------------------------------------------------------------------------#
# External dependencies
#------------------------------------------------------------------------------------#
if (ENABLE_MPI)
	find_package(MPI)
endif()

#------------------------------------------------------------------------------------#
# Customized build types
#------------------------------------------------------------------------------------#

# GNU Profiling
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

# Scalasca Profiling
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

# Set build type variable
if(NOT CMAKE_BUILD_TYPE)
	# Set default build type to Debug
	set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
		"Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel GNUProfiling ScalascaProfiling."
		FORCE)

	# Set the possible values of build type for the GUI
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
	"MinSizeRel" "RelWithDebInfo" "GNUProfiling" "ScalascaProfiling")
endif()

#------------------------------------------------------------------------------------#
# Compilation settings
#------------------------------------------------------------------------------------#
set (BITPIT_DEFINITIONS_PRIVATE "")
set (BITPIT_DEFINITIONS_PUBLIC "")

string(TOLOWER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_LOWER)

set (ENABLE_WARNINGS ${VERBOSE_MAKE})
if (VERBOSE_MAKE)
	set(CMAKE_VERBOSE_MAKEFILE ON)
endif()

if (ENABLE_MPI)
	list (APPEND BITPIT_DEFINITIONS_PUBLIC "BITPIT_ENABLE_MPI=1")

	if (MPI_C_FOUND)
		set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
	else()
		set(CMAKE_C_COMPILER mpicc)
	endif()

	if (MPI_CXX_FOUND)
		set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
	else()
		set(CMAKE_CXX_COMPILER mpic++)
	endif()
else ()
	list (APPEND BITPIT_DEFINITIONS_PUBLIC "BITPIT_ENABLE_MPI=0")
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g")
set(CMAKE_C_FLAGS_RELEASE "-O2")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

if (ENABLE_WARNINGS)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra") 
endif()

if (NOT ("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_LESS "2.8.12"))
	add_compile_options("-std=c++11")
else ()
	add_definitions("-std=c++11")
endif ()

# Define a preprocessor macro to recognize debug builds
IF(CMAKE_BUILD_TYPE_LOWER MATCHES "debug" OR CMAKE_BUILD_TYPE_LOWER MATCHES "debinfo")
	list (APPEND BITPIT_DEFINITIONS_PRIVATE "ENABLE_DEBUG=1")
else ()
	list (APPEND BITPIT_DEFINITIONS_PRIVATE "ENABLE_DEBUG=0")
endif ()

IF(NOT CMAKE_BUILD_TYPE_LOWER MATCHES "debug")
	list (APPEND BITPIT_DEFINITIONS_PRIVATE "NDEBUG")
endif ()

# Define an alias for building with scalasca
if (ENABLE_MPI)
	SET(C_FLAGS_INSTRUMENT   "-instrument mpicxx")
	SET(CXX_FLAGS_INSTRUMENT "-instrument mpic")
else ()
	SET(C_FLAGS_INSTRUMENT   "")
	SET(CXX_FLAGS_INSTRUMENT "")
endif ()

if (CMAKE_BUILD_TYPE_LOWER MATCHES "scalasca")
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

	file(INSTALL scalasca_cxx_compiler DESTINATION ${PROJECT_BINARY_DIR} PERMISSIONS OWNER_READ OWNER_EXECUTE )
	file(INSTALL scalasca_c_compiler   DESTINATION ${PROJECT_BINARY_DIR} PERMISSIONS OWNER_READ OWNER_EXECUTE )
	file(REMOVE  scalasca_cxx_compiler)
	file(REMOVE  scalasca_c_compiler)

	set(CMAKE_CXX_COMPILER "${PROJECT_BINARY_DIR}/scalasca_cxx_compiler")
	set(CMAKE_C_COMPILER   "${PROJECT_BINARY_DIR}/scalasca_c_compiler")
endif ()

# Add the definitions
set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS ${BITPIT_DEFINITIONS_PRIVATE})
set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS ${BITPIT_DEFINITIONS_PUBLIC})

#------------------------------------------------------------------------------------#
# Library version
#------------------------------------------------------------------------------------#
set (VERSION_HEADER_FILENAME "${PROJECT_SOURCE_DIR}/src/common/${CMAKE_PROJECT_NAME}_version.hpp")
string(TOUPPER "${CMAKE_PROJECT_NAME}_VERSION" VERSION_DEFINE_NAME)

file(READ "${VERSION_HEADER_FILENAME}" header)

string(REGEX REPLACE ".*#[ \t]*define[ \t]*${VERSION_DEFINE_NAME}[ \t]*\"([^\n]*)\".*" "\\1" match "${header}")
set(BITPIT_VERSION "${match}")

STRING(REGEX REPLACE "^([0-9]+)\\.[0-9]+\\.[0-9]+(-[0-9A-Za-z-]+)?" "\\1" match "${BITPIT_VERSION}")
set(BITPIT_MAJOR_VERSION "${match}")

STRING(REGEX REPLACE "^[0-9]+\\.([0-9])+\\.[0-9]+(-[0-9A-Za-z-]+)?" "\\1" match "${BITPIT_VERSION}")
set(BITPIT_MINOR_VERSION "${match}")

STRING(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+)(-[0-9A-Za-z-]+)?" "\\1" match "${BITPIT_VERSION}")
set(BITPIT_PATCH_VERSION "${match}")

STRING(REGEX MATCH "^[0-9]+\\.[0-9]+\\.[0-9]+-([0-9A-Za-z-]+)" match "${BITPIT_VERSION}")
if (NOT match STREQUAL "")
	STRING(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.[0-9]+-([0-9A-Za-z-]+)" "\\1" match "${BITPIT_VERSION}")
	set(BITPIT_TAG_VERSION "${match}")
else ()
	set(BITPIT_TAG_VERSION "")
endif ()

#------------------------------------------------------------------------------------#
# List of all avilable modules
#------------------------------------------------------------------------------------#
set(BITPIT_MODULE_LIST "")
list(APPEND BITPIT_MODULE_LIST "common")
list(APPEND BITPIT_MODULE_LIST "operators")
list(APPEND BITPIT_MODULE_LIST "containers")
list(APPEND BITPIT_MODULE_LIST "IO")
list(APPEND BITPIT_MODULE_LIST "LA")
list(APPEND BITPIT_MODULE_LIST "SA")
list(APPEND BITPIT_MODULE_LIST "CG")
list(APPEND BITPIT_MODULE_LIST "PABLO")
list(APPEND BITPIT_MODULE_LIST "patchkernel")
list(APPEND BITPIT_MODULE_LIST "surfunstructured")
list(APPEND BITPIT_MODULE_LIST "volcartesian")
list(APPEND BITPIT_MODULE_LIST "voloctree")
list(APPEND BITPIT_MODULE_LIST "volunstructured")
list(APPEND BITPIT_MODULE_LIST "RBF")
list(APPEND BITPIT_MODULE_LIST "levelset")

#------------------------------------------------------------------------------------#
# Modules dependecies
#------------------------------------------------------------------------------------#
set(OPERATORS_DEPS "common")
set(CONTAINERS_DEPS "common;operators")
set(IO_DEPS "common;operators")
set(LA_DEPS "common;operators")
set(SA_DEPS "common;operators")
set(CG_DEPS "common;operators")
set(PABLO_DEPS "IO")
set(PATCHKERNEL_DEPS "common;operators;containers;IO;LA;SA")
set(SURFUNSTRUCTURED_DEPS "common;patchkernel")
set(VOLCARTESIAN_DEPS "common;patchkernel")
set(VOLOCTREE_DEPS "common;PABLO;patchkernel")
set(VOLUNSTRUCTURED_DEPS "common;patchkernel")
set(RBF_DEPS "operators,SA")
set(LEVELSET_DEPS "common;SA;surfunstructured;voloctree,volcartesian")

#------------------------------------------------------------------------------------#
# Check of the dependencies
#------------------------------------------------------------------------------------#
foreach (MODULE_NAME IN LISTS BITPIT_MODULE_LIST)
	enableModule(${MODULE_NAME} TRUE FALSE)
endforeach ()

set(CHECK_DEPS 1)
while (${CHECK_DEPS} EQUAL 1)
	set(CHECK_DEPS 0)
	foreach (MODULE_NAME IN LISTS BITPIT_MODULE_LIST)
		isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)
		if (MODULE_ENABLED)
			string(TOUPPER ${MODULE_NAME} UPPER_MODULE_NAME)
			foreach (DEPENDENCY_NAME IN LISTS ${UPPER_MODULE_NAME}_DEPS)
				isModuleEnabled(${DEPENDENCY_NAME} DEPENDENCY_ENABLED)
				if (NOT DEPENDENCY_ENABLED)
					set(CHECK_DEPS 1)
					enableModule(${DEPENDENCY_NAME} TRUE TRUE)
					message("Module '${DEPENDENCY_NAME}' has been enabled because it's needed by module '${MODULE_NAME}'")
				endif ()
			endforeach ()
		endif ()
		unset(MODULE_ENABLE_FLAG)
	endforeach ()
endwhile()
unset(CHECK_DEPS)

# The modules 'common' and 'operators' has to be always built
getModuleEnableName("common" COMMON_MODULE_ENABLED)
unset(${COMMON_MODULE_ENABLED} CACHE)
set(${COMMON_MODULE_ENABLED} 1 CACHE INTERNAL "Request building common module" FORCE)

getModuleEnableName("operators" OPERATORS_MODULE_ENABLED)
unset(${OPERATORS_MODULE_ENABLED} CACHE)
set(${OPERATORS_MODULE_ENABLED} 1 CACHE INTERNAL "Request building operators module" FORCE)

#------------------------------------------------------------------------------------#
# External modules dependecies
#------------------------------------------------------------------------------------#

set(BITPIT_EXTERNAL_DEPENDENCIES "")

isModuleEnabled("CG" MODULE_CG_ENABLED)
isModuleEnabled("RBF" MODULE_RBF_ENABLED)
if (MODULE_CG_ENABLED OR MODULE_RBF_ENABLED)
	find_package(BLAS REQUIRED)
	find_package(LAPACK REQUIRED)

	foreach (LAPACK_LIBRARY_PATH IN LISTS LAPACK_LIBRARIES)
		get_filename_component(FILE_NAME "${LAPACK_LIBRARY_PATH}" NAME_WE)
		if ("${FILE_NAME}" STREQUAL "liblapack")
			get_filename_component(FILE_DIR "${LAPACK_LIBRARY_PATH}" PATH)
			get_filename_component(FILE_EXT "${LAPACK_LIBRARY_PATH}" EXT)

			file(TO_NATIVE_PATH "${FILE_DIR}/liblapacke${FILE_EXT}" LAPACKE_LIB_PATH)
			list (APPEND LAPACK_LIBRARIES "${LAPACKE_LIB_PATH}")
		endif ()
	endforeach ()

	list (APPEND BITPIT_EXTERNAL_DEPENDENCIES "${LAPACK_LIBRARIES}")
endif()

#------------------------------------------------------------------------------------#
# Subdirectories
#------------------------------------------------------------------------------------#

# Main library
add_subdirectory(src)

# CMake
add_subdirectory(cmake)

# Examples
add_subdirectory(examples)

# Tests
add_subdirectory(test)

# Docs
add_subdirectory(doc)

#------------------------------------------------------------------------------------#
# Tests
#------------------------------------------------------------------------------------#
enable_testing()

foreach (MODULE_NAME IN LISTS BITPIT_MODULE_LIST)
	isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)
	string(TOUPPER ${MODULE_NAME} UPPER_MODULE_NAME)
	if (MODULE_ENABLED AND DEFINED ${UPPER_MODULE_NAME}_TEST_ENTRIES)
		foreach(TEST_ENTRY IN LISTS ${UPPER_MODULE_NAME}_TEST_ENTRIES)
			getTestName("${TEST_ENTRY}" TEST_NAME)
			isTestParallel("${TEST_ENTRY}" IS_TEST_PARALLEL)
			if (${IS_TEST_PARALLEL})
				if (ENABLE_MPI)
					getTestNumProcs("${TEST_ENTRY}" "${MPIEXEC_MAX_NUMPROCS}" N_PROCS)
					add_test(NAME ${TEST_NAME} COMMAND ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${N_PROCS} ${MPIEXEC_POSTFLAGS} ${TEST_NAME} WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/test/${MODULE_NAME}")
				endif ()
			else ()
				add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME} WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/test/${MODULE_NAME}")
			endif ()
		endforeach ()
	endif ()
endforeach ()

add_custom_target(check DEPENDS tests COMMAND ${CMAKE_MAKE_PROGRAM} test)
