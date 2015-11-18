# CMake settings
cmake_minimum_required(VERSION 2.8)

# Add library to targets
include_directories("${PROJECT_SOURCE_DIR}/src/pablo")

file(GLOB SOURCE_FILES "*.cpp")
set(PABLO_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of Pablo Library" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp" "*.dat")
set(PABLO_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Headers of Pablo Library" FORCE)
unset(HEADER_FILES)

add_library(PABLO_TARGET_OBJECT OBJECT ${PABLO_SOURCES})

install(FILES ${PABLO_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
