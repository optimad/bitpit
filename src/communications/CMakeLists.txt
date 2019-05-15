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

set(COMMUNICATIONS_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of communications Library" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
set(COMMUNICATIONS_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Headers of communications Library" FORCE)
unset(HEADER_FILES)

if (NOT "${COMMUNICATIONS_SOURCES}" STREQUAL "")
    set(COMMUNICATIONS_TARGET_OBJECT "COMMUNICATIONS_TARGET_OBJECT")
    add_library(${COMMUNICATIONS_TARGET_OBJECT} OBJECT ${COMMUNICATIONS_SOURCES})
    set_lto_property(${COMMUNICATIONS_TARGET_OBJECT})
endif ()

install(FILES ${COMMUNICATIONS_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
