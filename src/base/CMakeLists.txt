# CMake settings
cmake_minimum_required(VERSION 2.8)

# Library name
set(BASE_LIBRARY "bitpit_base" CACHE INTERNAL "Base Library" FORCE)

# Library type
if (STATIC)
	set(LIBRARY_TYPE "STATIC")
else()
	set(LIBRARY_TYPE "SHARED")
endif()

# Add library to targets
file(GLOB SOURCE_FILES "*.cpp")
set(BASE_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of base Library" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
set(BASE_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Headers of base Library" FORCE)
unset(BASE_HEADERS)

add_library(${BASE_LIBRARY} ${LIBRARY_TYPE} ${BASE_SOURCES})
if (WITHOUT_MPI)
	set_target_properties(${BASE_LIBRARY} PROPERTIES DEBUG_POSTFIX "_D")
else()
	set_target_properties(${BASE_LIBRARY} PROPERTIES DEBUG_POSTFIX "_MPI_D")
	set_target_properties(${BASE_LIBRARY} PROPERTIES RELEASE_POSTFIX "_MPI")
endif()

install(TARGETS ${BASE_LIBRARY} DESTINATION lib)
install(FILES ${BASE_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
