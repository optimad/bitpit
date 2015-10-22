# Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

# List of tests
set(TESTS "")
list(APPEND TESTS "test_base_00001")
if (IO_ENABLE_VTK_WRAPPERS)
	list(APPEND TESTS "test_IO_00001")
endif()
list(APPEND TESTS "test_SA_00001")
list(APPEND TESTS "test_LA_00001")

set(TEST_LIST "${TESTS}" CACHE INTERNAL "List of tests" FORCE)

# External dependencies
if (IO_ENABLE_VTK_WRAPPERS)
	find_package(Boost ${IO_REQUIRED_BOOST_VERSION} REQUIRED)
endif()

# Rules to build the tests
include_directories("${PROJECT_SOURCE_DIR}/src/base")
include_directories("${PROJECT_SOURCE_DIR}/src/IO")
include_directories("${PROJECT_SOURCE_DIR}/src/SA")
include_directories("${PROJECT_SOURCE_DIR}/src/LA")
if (IO_ENABLE_VTK_WRAPPERS)
	include_directories("${Boost_INCLUDE_DIRS}")
endif()

foreach(TEST_NAME IN LISTS TEST_LIST)
	set(TEST_SOURCES "${TEST_NAME}.cpp")

	add_executable(${TEST_NAME} "${TEST_SOURCES}")
	target_link_libraries(${TEST_NAME} ${BITP_BASE_LIBRARY})
	if (IO_ENABLE_VTK_WRAPPERS)
		target_link_libraries(${TEST_NAME} "${Boost_LIBRARIES}")
	endif()
endforeach()

add_custom_target(tests DEPENDS ${TEST_LIST})
add_custom_target(clean-tests COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Data files
add_definitions(-DDATAFILE_PATH="${CMAKE_CURRENT_BINARY_DIR}/data/")

if (IO_ENABLE_VTK_WRAPPERS)
	add_custom_command(
		TARGET "test_IO_00001" PRE_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/selection.vtu" "${CMAKE_CURRENT_BINARY_DIR}/data/selection.vtu"
	)
endif()
