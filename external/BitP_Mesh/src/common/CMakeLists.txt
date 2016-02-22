# CMake settings
cmake_minimum_required(VERSION 2.8)

# Add common files to targets
file(GLOB SOURCE_FILES "*.cpp")
set(COMMON_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Common sources" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
set(COMMON_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Common headers" FORCE)
unset(HEADER_FILES)

add_library(COMMON_TARGET_OBJECT OBJECT ${COMMON_SOURCES})

install(FILES ${COMMON_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
