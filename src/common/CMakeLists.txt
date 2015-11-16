# CMake settings
cmake_minimum_required(VERSION 2.8)

# Add common files to targets
file(GLOB HEADER_FILES "*.hpp" "*.tpp")
set(COMMON_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Common headers" FORCE)
unset(HEADER_FILES)

install(FILES ${COMMON_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
