# CMake settings
cmake_minimum_required(VERSION 2.8)

# Add library to targets
file(GLOB SOURCE_FILES "*.cpp")
set(OPERATORS_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of the operators" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
set(OPERATORS_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Headers of the operators" FORCE)
unset(OPERATORS_HEADERS)

if (NOT "${OPERATORS_SOURCES}" STREQUAL "")
    add_library(OPERATORS_TARGET_OBJECT OBJECT ${OPERATORS_SOURCES})
endif ()

install(FILES ${OPERATORS_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
