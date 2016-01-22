# CMake settings
cmake_minimum_required(VERSION 2.8)

# Add library to targets
include_directories("${PROJECT_SOURCE_DIR}/src/operators")

file(GLOB SOURCE_FILES "*.cpp")
set(LA_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of LA Library" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
set(LA_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Headers of LA Library" FORCE)
unset(HEADER_FILES)

if (NOT "${LA_SOURCES}" STREQUAL "")
    add_library(LA_TARGET_OBJECT OBJECT ${LA_SOURCES})
endif ()

install(FILES ${LA_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
