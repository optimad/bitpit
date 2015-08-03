#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

#Name your project here
project("patchman")

#------------------------------------------------------------------------------------#
# Variables visible to the user
#------------------------------------------------------------------------------------#
set(ENABLE_PROFILING 0 CACHE BOOL "Set appropriate compiler flags to enable profiling fo GLORIA")
set(STATIC 1 CACHE BOOL "If set, the program is linked using static libraries")
set(BITPIT_LIB_PATH "/opt/bitpit" CACHE PATH "Directory with the BITPIT libraries")
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	set(CMAKE_INSTALL_PREFIX "/opt/bitpit/${CMAKE_PROJECT_NAME}" CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
endif()
#------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------#
# Internal variables
#------------------------------------------------------------------------------------#
set(WITHOUT_MPI "0" CACHE INTERNAL "If set, the program is compiled without MPI support")
#------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------#
# Compiler settings
#------------------------------------------------------------------------------------#
if (WITHOUT_MPI)
	add_definitions(-DDISABLE_MPI)
else()
	set(CMAKE_C_COMPILER mpicc)
	set(CMAKE_CXX_COMPILER mpic++)
endif()

if (ENABLE_PROFILING)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
	set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
endif (ENABLE_PROFILING)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -g -O0")

ADD_DEFINITIONS(-std=c++0x)
#------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------#
# Make settings
#------------------------------------------------------------------------------------#
if (CMAKE_BUILD_TYPE MATCHES Debug)
	set(CMAKE_VERBOSE_MAKEFILE ON)
endif (CMAKE_BUILD_TYPE MATCHES Debug)

#------------------------------------------------------------------------------------#
# Internal dependencies
#------------------------------------------------------------------------------------#
include_directories("include")
include_directories("src")
#------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------#
# External dependencies
#------------------------------------------------------------------------------------#
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake;${BITPIT_LIB_PATH}/cmake/Modules")

find_package(PABLO)
if (PABLO_FOUND)
	include_directories(${PABLO_INCLUDE_DIRS})
else()
	add_definitions(-DDISABLE_OCTREE)
endif()

set(VTK_DIR "/home/andrea/work/codes/libraries/VTK-6.2.0/VTK-6.2.0-build")
find_package(VTK 6.0 REQUIRED COMPONENTS vtkIOXML vtkIOLegacy NO_MODULE)
include(${VTK_USE_FILE})
#------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------#
# Subdirectories
#------------------------------------------------------------------------------------#

# Main program
add_subdirectory(src)

# Docs
add_subdirectory(doc)

# CMake
add_subdirectory(cmake)

# Tests
add_subdirectory(test)
#------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------#
# Tests
#------------------------------------------------------------------------------------#
enable_testing()

foreach(TEST_NAME IN LISTS TEST_LIST)
	add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()
#------------------------------------------------------------------------------------#

