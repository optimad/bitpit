# Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

# Name your project here
set(PABLO_LIBRARY "PABLO" CACHE INTERNAL "PABLO Library" FORCE)

# Library type
if (STATIC)
	set(LIBRARY_TYPE "STATIC")
else()
	set(LIBRARY_TYPE "SHARED")
endif()

# Add library to targets
file(GLOB SOURCE_FILES "*.cpp")
set(PABLO_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of PABLO Library" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp")
set(PABLO_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Headers of PABLO Library" FORCE)
unset(BASE_HEADERS)

add_library(${PABLO_LIBRARY} ${LIBRARY_TYPE} ${PABLO_SOURCES})
if (WITHOUT_MPI)
	set_target_properties(${PABLO_LIBRARY} PROPERTIES DEBUG_POSTFIX "_D")
else()
	set_target_properties(${PABLO_LIBRARY} PROPERTIES DEBUG_POSTFIX "_MPI_D")
	set_target_properties(${PABLO_LIBRARY} PROPERTIES RELEASE_POSTFIX "_MPI")
endif()

set_target_properties(${PABLO_LIBRARY} PROPERTIES VERSION "${PABLO_VERSION}"
                                                     SOVERSION  "${PABLO_MAJOR_VERSION}")

add_custom_target(library DEPENDS ${PABLO_LIBRARY})
add_custom_target(clean-library COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

install(TARGETS ${PABLO_LIBRARY} DESTINATION lib)
install(FILES ${PABLO_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})