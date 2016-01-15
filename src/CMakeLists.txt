#------------------------------------------------------------------------------------#
# Subdirectories
#------------------------------------------------------------------------------------#

# Common
add_subdirectory(common)

# Operators
add_subdirectory(operators)

# Containers
add_subdirectory(containers)

# IO
add_subdirectory(IO)

# SA
add_subdirectory(SA)

# LA
add_subdirectory(LA)

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

add_library(${BITP_BASE_LIBRARY} ${LIBRARY_TYPE}
	$<TARGET_OBJECTS:OPERATORS_TARGET_OBJECT>
	$<TARGET_OBJECTS:CONTAINERS_TARGET_OBJECT>
	$<TARGET_OBJECTS:IO_TARGET_OBJECT>
	$<TARGET_OBJECTS:SA_TARGET_OBJECT>
	$<TARGET_OBJECTS:LA_TARGET_OBJECT>
)

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
