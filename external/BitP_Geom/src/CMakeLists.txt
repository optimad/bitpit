#------------------------------------------------------------------------------------#
# Subdirectories
#------------------------------------------------------------------------------------#

# BASE
add_subdirectory(base)
include_directories("base")

# LEVELSET
add_subdirectory(levelset)

# DELAUNAY
#add_subdirectory(delaunay)

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
set(BITP_GEOM_LIBRARY ${CMAKE_PROJECT_NAME} CACHE INTERNAL "BitPit Geometry Library" FORCE)

add_library(${BITP_GEOM_LIBRARY} ${LIBRARY_TYPE}
    $<TARGET_OBJECTS:BASE_TARGET_OBJECT>
    $<TARGET_OBJECTS:LEVELSET_TARGET_OBJECT>
    #    $<TARGET_OBJECTS:DELAUNAY_TARGET_OBJECT>
)

if (ENABLE_MPI)
    set_target_properties(${BITP_GEOM_LIBRARY} PROPERTIES DEBUG_POSTFIX "_MPI_D")
    set_target_properties(${BITP_GEOM_LIBRARY} PROPERTIES RELEASE_POSTFIX "_MPI")

else()
    set_target_properties(${BITP_GEOM_LIBRARY} PROPERTIES DEBUG_POSTFIX "_D")
endif()

set_target_properties(${BITP_GEOM_LIBRARY} PROPERTIES VERSION "${BITP_GEOM_VERSION}"
    SOVERSION  "${BITP_GEOM_MAJOR_VERSION}")

install(TARGETS ${BITP_GEOM_LIBRARY} DESTINATION lib)

# Targets
add_custom_target(library DEPENDS ${GEOM_LIBRARY})
add_custom_target(clean-library COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
