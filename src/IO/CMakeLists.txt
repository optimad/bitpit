#---------------------------------------------------------------------------*\
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
#  as published by by the Free Software Foundation.
#
#  BitPit is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
#  for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
#
#*---------------------------------------------------------------------------*/

# CMake settings
cmake_minimum_required(VERSION 2.8)

# Features
set(IO_ENABLE_VTK_WRAPPERS 0 CACHE BOOL "Enable VTK wrappers (requires the library Boost)")

# Required libraries
if (IO_ENABLE_VTK_WRAPPERS)
	set(IO_REQUIRED_BOOST_VERSION "1.50.0" CACHE INTERNAL "Required version of the Boost library" FORCE)
	find_package(Boost ${IO_REQUIRED_BOOST_VERSION} REQUIRED)
endif()

# VTK Wrappers feature
set(IO_VTK_WRAPPERS_HEADER_HPP "VTKWrappers.hpp")
set(IO_VTK_WRAPPERS_HEADER_TPP "VTKWrappers.tpp")
set(IO_VTK_WRAPPERS_SOURCE_CPP "VTKWrappers.cpp")
set(IO_VTK_WRAPPERS_DEFINITION "IO_ENABLE_VTK_WRAPPERS")

# Add library to targets
include_directories("${PROJECT_SOURCE_DIR}/src/common")
include_directories("${PROJECT_SOURCE_DIR}/src/operators")
if (IO_ENABLE_VTK_WRAPPERS)
	include_directories("${Boost_INCLUDE_DIRS}")
endif()

file(GLOB SOURCE_FILES "*.cpp")
if (NOT IO_ENABLE_VTK_WRAPPERS)
	list(REMOVE_ITEM SOURCE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/${IO_VTK_WRAPPERS_SOURCE_CPP}")
endif()
set(IO_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of IO Library" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
if (NOT IO_ENABLE_VTK_WRAPPERS)
	list(REMOVE_ITEM HEADER_FILES "${CMAKE_CURRENT_SOURCE_DIR}/${IO_VTK_WRAPPERS_HEADER_HPP}")
	list(REMOVE_ITEM HEADER_FILES "${CMAKE_CURRENT_SOURCE_DIR}/${IO_VTK_WRAPPERS_HEADER_TPP}")
endif()
set(IO_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Headers of IO Library" FORCE)
unset(HEADER_FILES)

if (NOT "${IO_SOURCES}" STREQUAL "")
	add_library(IO_TARGET_OBJECT OBJECT ${IO_SOURCES})
endif ()

install(FILES ${IO_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})

# Definitions
set(IO_DEFINITIONS_DESCRIPTION "Pre-processor definitions needed by IO Library")
if (IO_ENABLE_VTK_WRAPPERS)
	add_definitions(-D${IO_VTK_WRAPPERS_DEFINITION})

	set(IO_DEFINITIONS "${IO_VTK_WRAPPERS_DEFINITION}" CACHE INTERNAL "${IO_DEFINITIONS_DESCRIPTION}" FORCE)
else ()
	set(IO_DEFINITIONS "" CACHE INTERNAL "${IO_DEFINITIONS_DESCRIPTION}" FORCE)
endif()
