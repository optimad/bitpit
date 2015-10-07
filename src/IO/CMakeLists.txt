# CMake settings
cmake_minimum_required(VERSION 2.8)

# Add library to targets
include_directories("${PROJECT_SOURCE_DIR}/src/base")

file(GLOB SOURCE_FILES "*.cpp")
set(IO_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of IO Library" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
set(IO_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Headers of IO Library" FORCE)
unset(HEADER_FILES)

add_library(IO_TARGET_OBJECT OBJECT ${IO_SOURCES})

install(FILES ${IO_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
