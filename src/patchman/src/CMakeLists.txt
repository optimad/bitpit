# Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

# Name your project here
set(PATCHMAN_LIBRARY "patchman" CACHE INTERNAL "PatchMan Library" FORCE)

# Library type
if (STATIC)
	set(LIBRARY_TYPE "STATIC")
else()
	set(LIBRARY_TYPE "SHARED")
endif()

# Add library to targets
file(GLOB SOURCE_FILES "*.cpp")
set(PATCHMAN_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of PatchMan Library" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp")
set(PATCHMAN_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Headers of PatchMan Library" FORCE)
unset(BASE_HEADERS)

add_library(${PATCHMAN_LIBRARY} ${LIBRARY_TYPE} ${PATCHMAN_SOURCES})
if (WITHOUT_MPI)
	set_target_properties(${PATCHMAN_LIBRARY} PROPERTIES DEBUG_POSTFIX "_D")
else()
	set_target_properties(${PATCHMAN_LIBRARY} PROPERTIES DEBUG_POSTFIX "_MPI_D")
	set_target_properties(${PATCHMAN_LIBRARY} PROPERTIES RELEASE_POSTFIX "_MPI")
endif()

set_target_properties(${PATCHMAN_LIBRARY} PROPERTIES VERSION "${PATCHMAN_VERSION}"
                                                     SOVERSION  "${PATCHMAN_MAJOR_VERSION}")

add_custom_target(library DEPENDS ${PATCHMAN_LIBRARY})
add_custom_target(clean-library COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

install(TARGETS ${PATCHMAN_LIBRARY} DESTINATION lib)
install(FILES ${PATCHMAN_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
