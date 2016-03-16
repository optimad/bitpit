# CMake settings
cmake_minimum_required(VERSION 2.8)

# Add library to targets
include_directories("${PROJECT_SOURCE_DIR}/src/base")

file(GLOB SOURCE_FILES "*.cpp")
set(LEVELSET_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of LEVELSET Library" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
set(LEVELSET_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Headers of LEVELSET Library" FORCE)
unset(HEADER_FILES)

add_library(LEVELSET_TARGET_OBJECT OBJECT ${LEVELSET_SOURCES})

install(FILES ${LEVELSET_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
