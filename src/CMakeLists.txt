# Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

# Name your project here
set(PATCHMAN_LIBRARY "patchman")
set (PATCHMAN_LIBRARY ${PATCHMAN_LIBRARY} PARENT_SCOPE)

# Library type
if (STATIC)
	set(LIBRARY_TYPE "STATIC")
else()
	set(LIBRARY_TYPE "SHARED")
endif()

# Library version
set (VERSION_HEADER_FILENAME "${CMAKE_CURRENT_SOURCE_DIR}/${PATCHMAN_LIBRARY}v.h")
set (VERSION_DEFINE_NAME "VERSION")

file(READ "${VERSION_HEADER_FILENAME}" header)

string(REGEX REPLACE ".*#[ \t]*define[ \t]*${VERSION_DEFINE_NAME}[ \t]*\"([^\n]*)\".*" "\\1" match "${header}")
set(PATCHMAN_VERSION "${match}")
set(PATCHMAN_VERSION "${PATCHMAN_VERSION}" PARENT_SCOPE)

STRING(REGEX REPLACE "^([0-9]+)\\.[0-9]+\\.[0-9]+" "\\1" match "${PATCHMAN_VERSION}")
set(PATCHMAN_MAJOR_VERSION "${match}")
set(PATCHMAN_MAJOR_VERSION "${PATCHMAN_MAJOR_VERSION}" PARENT_SCOPE)

# Add library to targets
file(GLOB PATCHMAN_SOURCES "*.cpp")
set (PATCHMAN_SOURCES ${PATCHMAN_SOURCES} PARENT_SCOPE)

file(GLOB PATCHMAN_HEADERS "*.hpp")
set (PATCHMAN_HEADERS ${PATCHMAN_HEADERS} PARENT_SCOPE)

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

install(TARGETS ${PATCHMAN_LIBRARY} DESTINATION lib)
install(FILES ${PATCHMAN_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
