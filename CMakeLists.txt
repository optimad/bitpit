#------------------------------------------------------------------------------------#
# CMake settings
#------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 0.99.0)

project("BitP_Base")

#------------------------------------------------------------------------------------#
# Variables visible to the user
#------------------------------------------------------------------------------------#
if (BITPIT_DIR)
	set(FIRST_BITPIT_DIR_INITIALIZATION 0)
else()
	set(FIRST_BITPIT_DIR_INITIALIZATION 1)
endif()

set(BITPIT_DIR "/opt/bitpit" CACHE PATH "The directory containing the BITPIT libraries")
if (${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT})
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

#------------------------------------------------------------------------------------#
# Internal variables
#------------------------------------------------------------------------------------#
set(WITHOUT_MPI "1" CACHE INTERNAL "If set, the program is compiled without MPI support")
set(STATIC 1 CACHE BOOL INTERNAL "If set, the program is linked using static libraries")

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

if (ENABLE_PROFILING)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
	set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
endif (ENABLE_PROFILING)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

if (ENABLE_WARNINGS)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -Wextra")
endif()

add_compile_options("-std=c++11")

#------------------------------------------------------------------------------------#
# External dependencies
#------------------------------------------------------------------------------------#
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake;${BITPIT_DIR}/cmake/Modules")

#------------------------------------------------------------------------------------#
# Subdirectories
#------------------------------------------------------------------------------------#

# Main library
add_subdirectory(src)

# opt_base
add_subdirectory(opt_base)

# opt_IO
add_subdirectory(opt_IO)

# opt_SA
add_subdirectory(opt_SA)

# opt_LA
add_subdirectory(opt_LA)

# Tests
add_subdirectory(test)

#------------------------------------------------------------------------------------#
# Tests
#------------------------------------------------------------------------------------#
enable_testing()

foreach(TEST_NAME IN LISTS TEST_LIST)
	add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()

add_custom_target(check DEPENDS tests COMMAND ${CMAKE_MAKE_PROGRAM} test)
