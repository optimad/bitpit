# CMake settings
cmake_minimum_required(VERSION 2.8)

# Library name
set(LA_LIBRARY "bitpit_LA" CACHE INTERNAL "LA Library" FORCE)

# Library type
if (STATIC)
	set(LIBRARY_TYPE "STATIC")
else()
	set(LIBRARY_TYPE "SHARED")
endif()

# Add library to targets
include_directories("${PROJECT_SOURCE_DIR}/src/base")

file(GLOB SOURCE_FILES "*.cpp")
set(LA_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of LA Library" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
set(LA_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Headers of LA Library" FORCE)
unset(HEADER_FILES)

add_library(${LA_LIBRARY} ${LIBRARY_TYPE} ${LA_SOURCES})
if (WITHOUT_MPI)
	set_target_properties(${LA_LIBRARY} PROPERTIES DEBUG_POSTFIX "_D")
else()
	set_target_properties(${LA_LIBRARY} PROPERTIES DEBUG_POSTFIX "_MPI_D")
	set_target_properties(${LA_LIBRARY} PROPERTIES RELEASE_POSTFIX "_MPI")
endif()

set_target_properties(${LA_LIBRARY} PROPERTIES VERSION "${BITP_BASE_VERSION}"
                                                 SOVERSION  "${BITP_BASE_MAJOR_VERSION}")

target_link_libraries(${LA_LIBRARY} ${BASE_LIBRARY})

install(TARGETS ${LA_LIBRARY} DESTINATION lib)
install(FILES ${LA_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
