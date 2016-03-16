# CMake settings
cmake_minimum_required(VERSION 2.8)

# Add library to targets
include_directories("${PROJECT_SOURCE_DIR}/src/base")

file(GLOB SOURCE_FILES "*.cpp")
set(BASE_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of BASE Library" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
set(BASE_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Headers of BASE Library" FORCE)
unset(HEADER_FILES)

add_library(BASE_TARGET_OBJECT OBJECT ${BASE_SOURCES})

install(FILES ${BASE_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
