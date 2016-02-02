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

# CMake settings
cmake_minimum_required(VERSION 2.8)

# Add library to targets
include_directories("${PROJECT_SOURCE_DIR}/src/common")
include_directories("${PROJECT_SOURCE_DIR}/src/operators")
include_directories("${PROJECT_SOURCE_DIR}/src/IO")

file(GLOB SOURCE_FILES "*.cpp")
set(UCARTMESH_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of ucartmesh module" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
set(UCARTMESH_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Headers of ucartmesh module" FORCE)
unset(HEADER_FILES)

if (NOT "${UCARTMESH_SOURCES}" STREQUAL "")
    add_library(UCARTMESH_TARGET_OBJECT OBJECT ${UCARTMESH_SOURCES})
endif ()

install(FILES ${UCARTMESH_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
