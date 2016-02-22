# CMake settings
cmake_minimum_required(VERSION 2.8)

# Add library to targets
include_directories("${PROJECT_SOURCE_DIR}/src/surftri")

file(GLOB SOURCE_FILES "*.cpp")
set(SURFTRI_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of SurfTri Library" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
set(SURFTRI_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Headers of SurfTri Library" FORCE)
unset(HEADER_FILES)

add_library(SURFTRI_TARGET_OBJECT OBJECT ${SURFTRI_SOURCES})

install(FILES ${SURFTRI_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
