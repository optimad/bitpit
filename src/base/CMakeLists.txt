# CMake settings
cmake_minimum_required(VERSION 2.8)

# Add library to targets
file(GLOB SOURCE_FILES "*.cpp")
set(BASE_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of base Library" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
set(BASE_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Headers of base Library" FORCE)
unset(BASE_HEADERS)

add_library(BASE_TARGET_OBJECT OBJECT ${BASE_SOURCES})

install(FILES ${BASE_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
