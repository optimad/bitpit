# CMake settings
cmake_minimum_required(VERSION 2.8)

# Add library to targets
include_directories("${PROJECT_SOURCE_DIR}/src/base")

file(GLOB SOURCE_FILES "*.cpp")
set(DELAUNAY_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of DELAUNAY Library" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
set(DELAUNAY_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Headers of DELAUNAY Library" FORCE)
unset(HEADER_FILES)

add_library(DELAUNAY_TARGET_OBJECT OBJECT ${DELAUNAY_SOURCES})

install(FILES ${DELAUNAY_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
