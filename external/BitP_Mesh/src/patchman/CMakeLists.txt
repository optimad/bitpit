# CMake settings
cmake_minimum_required(VERSION 2.8)

# Add library to targets
include_directories("${PROJECT_SOURCE_DIR}/src/common")

file(GLOB SOURCE_FILES "*.cpp")
set(PATCHMAN_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of PatchMan Library" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
set(PATCHMAN_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Headers of PatchMan Library" FORCE)
unset(HEADER_FILES)

add_library(PATCHMAN_TARGET_OBJECT OBJECT ${PATCHMAN_SOURCES})

install(FILES ${PATCHMAN_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
