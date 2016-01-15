# Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

# List of tests
set(TESTS "")
list(APPEND TESTS "test_LA_00001")

set(LA_TEST_LIST "${TESTS}" CACHE INTERNAL "List of tests for the LA module" FORCE)

# Rules to build the tests
include_directories("${PROJECT_SOURCE_DIR}/src/operators")
include_directories("${PROJECT_SOURCE_DIR}/src/containers")
include_directories("${PROJECT_SOURCE_DIR}/src/LA")

foreach(TEST_NAME IN LISTS LA_TEST_LIST)
	set(TEST_SOURCES "${TEST_NAME}.cpp")

	add_executable(${TEST_NAME} "${TEST_SOURCES}")
	target_link_libraries(${TEST_NAME} ${BITP_BASE_LIBRARY})
endforeach()

add_custom_target(tests-LA DEPENDS ${LA_TEST_LIST})
add_custom_target(clean-tests-LA COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
