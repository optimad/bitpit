# Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

# List of tests
set(TESTS "")
list(APPEND TESTS "test0")
list(APPEND TESTS "test1")
list(APPEND TESTS "test104")
if (NOT WITHOUT_MPI)
	set(PTESTS "")
	list(APPEND PTESTS "test12")
endif()

set(TEST_LIST "${TESTS}" CACHE INTERNAL "List of tests" FORCE)
set(PTEST_LIST "${PTESTS}" CACHE INTERNAL "List of parallel tests" FORCE)

# External dependencies
find_package(Boost ${IO_REQUIRED_BOOST_VERSION} REQUIRED)

# Rules to build the tests
include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${Boost_INCLUDE_DIRS}")

foreach(TEST_NAME IN LISTS TEST_LIST)
	set(TEST_SOURCES "${TEST_NAME}.cpp")
	add_executable(${TEST_NAME} "${TEST_SOURCES}")
	target_link_libraries(${TEST_NAME} ${PABLO_LIBRARY})
	target_link_libraries(${TEST_NAME} "${Boost_LIBRARIES}")
endforeach()
if (NOT WITHOUT_MPI)
	foreach(PTEST_NAME IN LISTS PTEST_LIST)
		set(PTEST_SOURCES "${PTEST_NAME}.cpp")
		add_executable(${PTEST_NAME} "${PTEST_SOURCES}")
		target_link_libraries(${PTEST_NAME} ${PABLO_LIBRARY})
		target_link_libraries(${PTEST_NAME} "${Boost_LIBRARIES}")
	endforeach()
endif()

add_custom_target(tests DEPENDS ${TEST_LIST} ${PTEST_LIST})

add_custom_target(clean-tests COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})




