#---------------------------------------------------------------------------
#
#  bitpit
#
#  Copyright (C) 2015-2019 OPTIMAD engineering Srl
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
#---------------------------------------------------------------------------

# CMake settings
cmake_minimum_required(VERSION 2.8)
initialize_lto_property()

# Name of the current module
get_filename_component(MODULE_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)

# Add library to targets
addModuleIncludeDirectories(${MODULE_NAME})

file(GLOB SOURCE_FILES "*.cpp")

set(IO_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of IO Library" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
set(IO_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Headers of IO Library" FORCE)
unset(HEADER_FILES)

if (NOT "${IO_SOURCES}" STREQUAL "")
    set(IO_TARGET_OBJECT "IO_TARGET_OBJECT")
    add_library(${IO_TARGET_OBJECT} OBJECT ${IO_SOURCES})
    set_lto_property(${IO_TARGET_OBJECT})
endif ()

install(FILES ${IO_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})

# Definitions
set(IO_DEFINITIONS_DESCRIPTION "Pre-processor definitions needed by IO Library")
set(IO_DEFINITIONS "" CACHE INTERNAL "${IO_DEFINITIONS_DESCRIPTION}" FORCE)
