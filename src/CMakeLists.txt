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
#----------------------------------------------------------------------------

#------------------------------------------------------------------------------------#
# Subdirectories
#------------------------------------------------------------------------------------#

if (DEFINED MODULE_SOURCES)
	unset(${MODULE_SOURCES})
	unset(${MODULE_DEFINITIONS})
endif ()

# Common
add_subdirectory(common)
set(MODULE_SOURCES ${MODULE_SOURCES} $<TARGET_OBJECTS:COMMON_TARGET_OBJECT>)

# Modules
foreach(MODULE_NAME IN LISTS BITPIT_MODULE_LIST)
	isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)
	if (MODULE_ENABLED)
		string(TOUPPER ${MODULE_NAME} UPPERCASE_MODULE_NAME)
		add_subdirectory(${MODULE_NAME})
		if (TARGET ${UPPERCASE_MODULE_NAME}_TARGET_OBJECT)
			set(MODULE_SOURCES ${MODULE_SOURCES} $<TARGET_OBJECTS:${UPPERCASE_MODULE_NAME}_TARGET_OBJECT>)
		endif ()
		if (DEFINED ${UPPERCASE_MODULE_NAME}_DEFINITIONS)
			set(MODULE_DEFINITIONS ${MODULE_DEFINITIONS} ${${UPPERCASE_MODULE_NAME}_DEFINITIONS})
		endif ()
	endif ()
endforeach()

#------------------------------------------------------------------------------------#
# Library
#------------------------------------------------------------------------------------#

# Library type
if (STATIC)
	set(LIBRARY_TYPE "STATIC")
else()
	set(LIBRARY_TYPE "SHARED")
endif()

# Library
set(BITPIT_LIBRARY ${CMAKE_PROJECT_NAME} CACHE INTERNAL "bitpit Base Library" FORCE)

add_library(${BITPIT_LIBRARY} ${LIBRARY_TYPE} ${MODULE_SOURCES})

if (ENABLE_MPI)
	set_target_properties(${BITPIT_LIBRARY} PROPERTIES DEBUG_POSTFIX "_MPI_D")
	set_target_properties(${BITPIT_LIBRARY} PROPERTIES RELWITHDEBINFO_POSTFIX "_MPI_D")
	set_target_properties(${BITPIT_LIBRARY} PROPERTIES RELEASE_POSTFIX "_MPI")
else()
	set_target_properties(${BITPIT_LIBRARY} PROPERTIES DEBUG_POSTFIX "_D")
	set_target_properties(${BITPIT_LIBRARY} PROPERTIES RELWITHDEBINFO_POSTFIX "_D")
endif()

set_target_properties(${BITPIT_LIBRARY} PROPERTIES VERSION "${BITPIT_VERSION}"
                                                 SOVERSION  "${BITPIT_MAJOR_VERSION}")

install(TARGETS ${BITPIT_LIBRARY} DESTINATION lib)

# Definitions
set(BITPIT_DEFINITIONS ${MODULE_DEFINITIONS} CACHE INTERNAL "Pre-processor definitions needed by the Base Library" FORCE)

# Targets
add_custom_target(library DEPENDS ${BASE_LIBRARY})
add_custom_target(clean-library COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
