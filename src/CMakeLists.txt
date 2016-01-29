#------------------------------------------------------------------------------------#
# Subdirectories
#------------------------------------------------------------------------------------#

# Common
add_subdirectory(common)

# PABLO
add_subdirectory(pablo)
include_directories("pablo")

if (NOT ONLY_PABLO)
# UCARTMESH
add_subdirectory(ucartmesh)

# SURFTRI
add_subdirectory(surftri)
include_directories("surftri")

# VOLFTRI
add_subdirectory(voltri)

# PATCHMAN
add_subdirectory(patchman)
endif()

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
set(BITP_MESH_LIBRARY ${CMAKE_PROJECT_NAME} CACHE INTERNAL "bitpit Mesh Library" FORCE)

if (ONLY_PABLO)
    add_library(${BITP_MESH_LIBRARY} ${LIBRARY_TYPE}
        $<TARGET_OBJECTS:PABLO_TARGET_OBJECT>
    )
else()
    add_library(${BITP_MESH_LIBRARY} ${LIBRARY_TYPE}
        $<TARGET_OBJECTS:COMMON_TARGET_OBJECT>
        $<TARGET_OBJECTS:PABLO_TARGET_OBJECT>
        $<TARGET_OBJECTS:UCARTMESH_TARGET_OBJECT>
        $<TARGET_OBJECTS:SURFTRI_TARGET_OBJECT>
        $<TARGET_OBJECTS:VOLTRI_TARGET_OBJECT>
        $<TARGET_OBJECTS:PATCHMAN_TARGET_OBJECT>
    )
endif()

if (NOT ENABLE_MPI)
	set_target_properties(${BITP_MESH_LIBRARY} PROPERTIES DEBUG_POSTFIX "_D")
else()
	set_target_properties(${BITP_MESH_LIBRARY} PROPERTIES DEBUG_POSTFIX "_MPI_D")
	set_target_properties(${BITP_MESH_LIBRARY} PROPERTIES RELEASE_POSTFIX "_MPI")
endif()

set_target_properties(${BITP_MESH_LIBRARY} PROPERTIES VERSION "${BITP_MESH_VERSION}"
                                                 SOVERSION  "${BITP_MESH_MAJOR_VERSION}")

install(TARGETS ${BITP_MESH_LIBRARY} DESTINATION lib)

# Targets
add_custom_target(library DEPENDS ${MESH_LIBRARY})
add_custom_target(clean-library COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
