# Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

#------------------------------------------------------------------------------------#
# Subdirectories
#------------------------------------------------------------------------------------#

set(TEST_LIST_DESCRIPTION "List of tests")

# Operators
add_subdirectory(operators)
set(TEST_LIST "${OPERATORS_TEST_LIST}" CACHE INTERNAL ${TEST_LIST_DESCRIPTION} FORCE)

# IO
add_subdirectory(IO)
set(TEST_LIST "${TEST_LIST};${IO_TEST_LIST}" CACHE INTERNAL ${TEST_LIST_DESCRIPTION} FORCE)

# SA
add_subdirectory(SA)
set(TEST_LIST "${TEST_LIST};${SA_TEST_LIST}" CACHE INTERNAL ${TEST_LIST_DESCRIPTION} FORCE)

# LA
add_subdirectory(LA)
set(TEST_LIST "${TEST_LIST};${LA_TEST_LIST}" CACHE INTERNAL ${TEST_LIST_DESCRIPTION} FORCE)

#------------------------------------------------------------------------------------#
# Targets
#------------------------------------------------------------------------------------#

add_custom_target(tests DEPENDS ${TEST_LIST})
add_custom_target(clean-tests COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
