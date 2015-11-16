#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

# List of tests
set(TESTS "")
list(APPEND TESTS "ucartmesh")
list(APPEND TESTS "surftri")
list(APPEND TESTS "voltri")
list(APPEND TESTS "pablo")
list(APPEND TESTS "patchman")

set(TEST_LIST "${TESTS}" CACHE INTERNAL "List of tests" FORCE)

#Rules to build the tests
foreach(TEST_NAME IN LISTS TEST_LIST)
	set(TEST_SOURCES "${TEST_NAME}.cpp")

	add_executable(${TEST_NAME} "${TEST_SOURCES}")
    target_link_libraries(${TEST_NAME} ${BITP_MESH_LIBRARY})
	target_link_libraries(${TEST_NAME} ${VTK_LIBRARIES})

endforeach()

add_custom_target(tests DEPENDS ${TEST_LIST})
add_custom_target(clean-tests COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
