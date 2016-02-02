# CMake settings
cmake_minimum_required(VERSION 2.8)

# Add library to targets
file(GLOB SOURCE_FILES "*.cpp")
set(PABLO_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of PABLO module" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
set(PABLO_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Headers of PABLO module" FORCE)
unset(HEADER_FILES)

if (NOT "${PABLO_SOURCES}" STREQUAL "")
    add_library(PABLO_TARGET_OBJECT OBJECT ${PABLO_SOURCES})
endif ()

install(FILES ${PABLO_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
