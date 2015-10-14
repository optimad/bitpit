# CMake settings
cmake_minimum_required(VERSION 2.8)

# Features
set(IO_ENABLE_VTK_WRAPPERS 0 CACHE BOOL "Enable VTK wrappers (requires the library Boost)")

# Required libraries
if (IO_ENABLE_VTK_WRAPPERS)
	find_package(Boost 1.50.0 REQUIRED)
endif()

# VTK Wrappers feature
set(IO_VTK_WRAPPERS_HEADER_HPP "Class_VTK_Wrappers.hpp" CACHE INTERNAL "Header file needed by the VTK Wrapper feature" FORCE)
set(IO_VTK_WRAPPERS_HEADER_TPP "Class_VTK_Wrappers.tpp" CACHE INTERNAL "Template header file needed by the VTK Wrapper feature" FORCE)
set(IO_VTK_WRAPPERS_SOURCE_CPP "Class_VTK_Wrappers.cpp" CACHE INTERNAL "Source file needed by the VTK Wrapper feature" FORCE)
set(IO_VTK_WRAPPERS_DEFINITION "-DIO_ENABLE_VTK_WRAPPERS" CACHE INTERNAL "Preprocessor defnition needed by the VTK Wrapper feature" FORCE)

# Add library to targets
include_directories("${PROJECT_SOURCE_DIR}/src/base")

file(GLOB SOURCE_FILES "*.cpp")
if (NOT IO_ENABLE_VTK_WRAPPERS)
	list(REMOVE_ITEM SOURCE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/${IO_VTK_WRAPPERS_SOURCE_CPP}")
endif()
set(IO_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of IO Library" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
if (NOT IO_ENABLE_VTK_WRAPPERS)
	list(REMOVE_ITEM HEADER_FILES "${CMAKE_CURRENT_SOURCE_DIR}/${IO_VTK_WRAPPERS_HEADER_HPP}")
	list(REMOVE_ITEM HEADER_FILES "${CMAKE_CURRENT_SOURCE_DIR}/${IO_VTK_WRAPPERS_HEADER_TPP}")
endif()
set(IO_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Headers of IO Library" FORCE)
unset(HEADER_FILES)

add_library(IO_TARGET_OBJECT OBJECT ${IO_SOURCES})

install(FILES ${IO_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})

if (IO_ENABLE_VTK_WRAPPERS)
	add_definitions(${IO_VTK_WRAPPERS_DEFINITION})
endif()
