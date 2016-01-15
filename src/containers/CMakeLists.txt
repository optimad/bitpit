# CMake settings
cmake_minimum_required(VERSION 2.8)

# Add library to targets
include_directories("${PROJECT_SOURCE_DIR}/src/containers")

file(GLOB SOURCE_FILES "*.cpp")
set(CONTAINERS_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources for BitPit custom containers" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
set(CONTAINERS_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Headers for BitPit custom containers" FORCE)
unset(HEADER_FILES)

add_library(CONTAINERS_TARGET_OBJECT OBJECT ${CONTAINERS_SOURCES})

install(FILES ${CONTAINERS_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
