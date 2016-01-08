#------------------------------------------------------------------------------------#
# Subdirectories
#------------------------------------------------------------------------------------#

if (DEFINED BITPIT_SOURCES)
	unset(${BITPIT_SOURCES})
endif ()

# Common
add_subdirectory(common)

# Modules
foreach(MODULE_NAME IN LISTS BITPIT_MODULE_LIST)
	string(TOUPPER ${MODULE_NAME} UPPERCASE_MODULE_NAME)

	add_subdirectory(${MODULE_NAME})
	set(BITPIT_SOURCES ${BITPIT_SOURCES} $<TARGET_OBJECTS:${UPPERCASE_MODULE_NAME}_TARGET_OBJECT>)
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
set(BITP_BASE_LIBRARY ${CMAKE_PROJECT_NAME} CACHE INTERNAL "BitPit Base Library" FORCE)

add_library(${BITP_BASE_LIBRARY} ${LIBRARY_TYPE} ${BITPIT_SOURCES})

if (WITHOUT_MPI)
	set_target_properties(${BITP_BASE_LIBRARY} PROPERTIES DEBUG_POSTFIX "_D")
else()
	set_target_properties(${BITP_BASE_LIBRARY} PROPERTIES DEBUG_POSTFIX "_MPI_D")
	set_target_properties(${BITP_BASE_LIBRARY} PROPERTIES RELEASE_POSTFIX "_MPI")
endif()

set_target_properties(${BITP_BASE_LIBRARY} PROPERTIES VERSION "${BITP_BASE_VERSION}"
                                                 SOVERSION  "${BITP_BASE_MAJOR_VERSION}")

install(TARGETS ${BITP_BASE_LIBRARY} DESTINATION lib)

# Targets
add_custom_target(library DEPENDS ${BASE_LIBRARY})
add_custom_target(clean-library COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
