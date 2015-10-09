#------------------------------------------------------------------------------------#
# CMake settings
#------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 2.8.12)

project("patchman")

#------------------------------------------------------------------------------------#
# Variables visible to the user
#------------------------------------------------------------------------------------#
if (BITPIT_DIR)
	set(FIRST_BITPIT_DIR_INITIALIZATION 0)
else()
	set(FIRST_BITPIT_DIR_INITIALIZATION 1)
endif()

set(BITPIT_DIR "/opt/bitpit" CACHE PATH "The directory containing the BITPIT libraries")
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	set(CMAKE_INSTALL_PREFIX "${BITPIT_DIR}/${CMAKE_PROJECT_NAME}" CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
endif()

if (${FIRST_BITPIT_DIR_INITIALIZATION})
	return()
else()
	if (NOT EXISTS ${BITPIT_DIR})
		set(errorMessage "REQUIRED DIRECTORY NOT FOUND\n We could not find the BITPIT directory. You need to adjust CMake configuration in order to continue building ${CMAKE_PROJECT_NAME}.")
		message(FATAL_ERROR "${errorMessage}")
	endif()
endif()

set(ENABLE_PROFILING 0 CACHE BOOL "Set appropriate compiler flags to enable profiling")
set(VERBOSE_MAKE 0 CACHE BOOL "Set appropriate compiler and cmake flags to enable verbose output from compilation")
set(STATIC 1 CACHE BOOL "If set, the program is linked using static libraries")
set(VTK_DIR "/opt/VTK/lib/cmake/vtk" CACHE PATH "The directory containing a CMake configuration file for VTK")

#------------------------------------------------------------------------------------#
# Internal variables
#------------------------------------------------------------------------------------#
set(WITHOUT_MPI "0" CACHE INTERNAL "If set, the program is compiled without MPI support")

#------------------------------------------------------------------------------------#
# Compilation settings
#------------------------------------------------------------------------------------#
set (ENABLE_WARNINGS ${VERBOSE_MAKE})
if (VERBOSE_MAKE)
	set(CMAKE_VERBOSE_MAKEFILE ON)
endif()

if (WITHOUT_MPI)
	add_definitions(-DDISABLE_MPI)
else()
	set(CMAKE_C_COMPILER mpicc)
	set(CMAKE_CXX_COMPILER mpic++)
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

add_compile_options("-std=c++11")

#------------------------------------------------------------------------------------#
# Internal dependencies
#------------------------------------------------------------------------------------#
include_directories("include")
include_directories("src")

#------------------------------------------------------------------------------------#
# External dependencies
#------------------------------------------------------------------------------------#
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake;${BITPIT_DIR}/cmake/Modules")

find_package(PABLO)
if (PABLO_FOUND)
	include_directories(${PABLO_INCLUDE_DIRS})
else()
	add_definitions(-DDISABLE_OCTREE)
endif()

find_package(VTK 6.0 REQUIRED COMPONENTS vtkIOXML vtkIOLegacy NO_MODULE)
include(${VTK_USE_FILE})

#------------------------------------------------------------------------------------#
# Library version
#------------------------------------------------------------------------------------#
set (VERSION_HEADER_FILENAME "${CMAKE_CURRENT_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}_version.hpp")
string(TOUPPER "${CMAKE_PROJECT_NAME}_VERSION" VERSION_DEFINE_NAME)

file(READ "${VERSION_HEADER_FILENAME}" header)

string(REGEX REPLACE ".*#[ \t]*define[ \t]*${VERSION_DEFINE_NAME}[ \t]*\"([^\n]*)\".*" "\\1" match "${header}")
set(PATCHMAN_VERSION "${match}")

STRING(REGEX REPLACE "^([0-9]+)\\.[0-9]+\\.[0-9]+" "\\1" match "${PATCHMAN_VERSION}")
set(PATCHMAN_MAJOR_VERSION "${match}")

#------------------------------------------------------------------------------------#
# Subdirectories
#------------------------------------------------------------------------------------#

# Main program
add_subdirectory(src)

# CMake
add_subdirectory(cmake)

# Tests
add_subdirectory(test)

# Examples
add_subdirectory(examples)

# Docs
add_subdirectory(doc)

#------------------------------------------------------------------------------------#
# Tests
#------------------------------------------------------------------------------------#
enable_testing()

foreach(TEST_NAME IN LISTS TEST_LIST)
	add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()

add_custom_target(check DEPENDS tests COMMAND ${CMAKE_MAKE_PROGRAM} test)
