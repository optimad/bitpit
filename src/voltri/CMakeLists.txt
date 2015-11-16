# CMake settings
cmake_minimum_required(VERSION 2.8)

# Add library to targets
include_directories("${PROJECT_SOURCE_DIR}/src/voltri")

file(GLOB SOURCE_FILES "*.cpp")
set(VOLTRI_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of VolTri Library" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
set(VOLTRI_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Headers of VolTri Library" FORCE)
unset(HEADER_FILES)

add_library(VOLTRI_TARGET_OBJECT OBJECT ${VOLTRI_SOURCES})

install(FILES ${VOLTRI_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
