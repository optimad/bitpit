# CMake settings
cmake_minimum_required(VERSION 2.8)

# Add library to targets
include_directories("${PROJECT_SOURCE_DIR}/src/base")

file(GLOB SOURCE_FILES "*.cpp")
set(SA_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of SA Library" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
set(SA_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Headers of SA Library" FORCE)
unset(HEADER_FILES)

add_library(SA_TARGET_OBJECT OBJECT ${SA_SOURCES})

install(FILES ${SA_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
