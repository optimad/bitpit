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
set(ENABLE_PROFILING 0 CACHE BOOL "Set appropriate compiler flags to enable profiling")
set(ENABLE_MPI 0 CACHE BOOL "If set, the program is compiled with MPI support")
set(VERBOSE_MAKE 0 CACHE BOOL "Set appropriate compiler and cmake flags to enable verbose output from compilation")

#------------------------------------------------------------------------------------#
# Internal variables
#------------------------------------------------------------------------------------#
set(STATIC 1 CACHE BOOL INTERNAL "If set, the program is linked using static libraries")

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

#------------------------------------------------------------------------------------#
# External dependencies
#------------------------------------------------------------------------------------#
if (ENABLE_MPI)
	find_package(MPI)
endif()

#------------------------------------------------------------------------------------#
# Compilation settings
#------------------------------------------------------------------------------------#
set (ENABLE_WARNINGS ${VERBOSE_MAKE})
if (VERBOSE_MAKE)
	set(CMAKE_VERBOSE_MAKEFILE ON)
endif()

if (ENABLE_MPI)
	add_definitions(-DENABLE_MPI=1)

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

if (ENABLE_PROFILING)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
	set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
endif (ENABLE_PROFILING)

if (NOT ("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_LESS "2.8.12"))
	add_compile_options("-std=c++11")
else ()
	add_definitions("-std=c++11")
endif ()

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
list(APPEND BITPIT_MODULE_LIST "operators")
list(APPEND BITPIT_MODULE_LIST "containers")
list(APPEND BITPIT_MODULE_LIST "IO")
list(APPEND BITPIT_MODULE_LIST "LA")
list(APPEND BITPIT_MODULE_LIST "SA")
list(APPEND BITPIT_MODULE_LIST "PABLO")
list(APPEND BITPIT_MODULE_LIST "patch")
list(APPEND BITPIT_MODULE_LIST "octreepatch")
list(APPEND BITPIT_MODULE_LIST "ucartmesh")

#------------------------------------------------------------------------------------#
# Modules dependecies
#------------------------------------------------------------------------------------#
set(CONTAINERS_DEPS "operators")
set(IO_DEPS "operators")
set(LA_DEPS "operators")
set(SA_DEPS "operators")
set(UCARTMESH_DEPS "IO")
set(PATCH_DEPS "operators;containers;IO;LA")
set(OCTREEPATCH_DEPS "PABLO;patch")

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

# Module 'operators' module has to be always build
getModuleEnableName("operators" BASE_ENABLE_VARIABLE)
unset(${BASE_ENABLE_VARIABLE} CACHE)
set(${BASE_ENABLE_VARIABLE} 1 CACHE INTERNAL "Request building operators module" FORCE)

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

foreach(MODULE_NAME IN LISTS BITPIT_MODULE_LIST)
	isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)
	if (MODULE_ENABLED)
		string(TOUPPER ${MODULE_NAME} UPPER_MODULE_NAME)

		# Serial tests
		if (DEFINED ${UPPER_MODULE_NAME}_TEST_LIST_SERIAL)
			foreach(TEST_NAME IN LISTS ${UPPER_MODULE_NAME}_TEST_LIST_SERIAL)
				add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME} WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/test/${MODULE_NAME}")
			endforeach()
		endif()

		# Parallel tests
		if (ENABLE_MPI)
			if (DEFINED ${UPPER_MODULE_NAME}_TEST_LIST_PARALLEL)
				foreach(TEST_NAME IN LISTS ${UPPER_MODULE_NAME}_TEST_LIST_PARALLEL)
					add_test(NAME ${TEST_NAME} COMMAND ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_POSTFLAGS} ${TEST_NAME} WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/test/${MODULE_NAME}")
				endforeach()
			endif()
		endif ()
	endif ()
endforeach()

add_custom_target(check DEPENDS tests COMMAND ${CMAKE_MAKE_PROGRAM} test)
