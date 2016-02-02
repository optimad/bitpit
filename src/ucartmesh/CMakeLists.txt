# CMake settings
cmake_minimum_required(VERSION 2.8)

# Add library to targets
include_directories("${PROJECT_SOURCE_DIR}/src/common")
include_directories("${PROJECT_SOURCE_DIR}/src/operators")
include_directories("${PROJECT_SOURCE_DIR}/src/IO")

file(GLOB SOURCE_FILES "*.cpp")
set(UCARTMESH_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of ucartmesh module" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
set(UCARTMESH_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Headers of ucartmesh module" FORCE)
unset(HEADER_FILES)

if (NOT "${UCARTMESH_SOURCES}" STREQUAL "")
    add_library(UCARTMESH_TARGET_OBJECT OBJECT ${UCARTMESH_SOURCES})
endif ()

install(FILES ${UCARTMESH_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
