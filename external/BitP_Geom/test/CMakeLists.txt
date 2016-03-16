#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

# List of tests
set(SERIAL_TESTS "")

list(APPEND SERIAL_TESTS "CGBase_001")
list(APPEND SERIAL_TESTS "CGBase_002")

list(APPEND SERIAL_TESTS "levelset_001")
list(APPEND SERIAL_TESTS "levelset_002")
#list(APPEND SERIAL_TESTS "levelset_003")
list(APPEND SERIAL_TESTS "levelset_004")

if (ENABLE_MPI)
    set(PARALLEL_TEST "")
    #    list(APPEND PARALLEL_TESTS "parallel_pablo_001")
endif()

set(TEST_LIST "${SERIAL_TESTS}" CACHE INTERNAL "List of serial tests" FORCE)
set(PARALLELTEST_LIST "${PARALLEL_TESTS}" CACHE INTERNAL "List of parallel tests" FORCE)

include_directories("${PROJECT_SOURCE_DIR}/src/common")
include_directories("${PROJECT_SOURCE_DIR}/src/base")
include_directories("${PROJECT_SOURCE_DIR}/src/levelset")
include_directories("${PROJECT_SOURCE_DIR}/src/delaunay")

include_directories(${BITP_GEOM_INCLUDE_DIRS})

#Rules to build the tests
foreach(TEST_NAME IN LISTS TEST_LIST PARALLELTEST_LIST)
    set(TEST_SOURCES "${TEST_NAME}.cpp")

    add_executable(${TEST_NAME} "${TEST_SOURCES}")
    target_link_libraries(${TEST_NAME} ${BITP_GEOM_LIBRARY})
    target_link_libraries(${TEST_NAME} ${BITPIT_LIBRARY})
    target_link_libraries(${TEST_NAME} "lapack")
    target_link_libraries(${TEST_NAME} "lapacke")

endforeach()

add_custom_target(tests DEPENDS ${TEST_LIST} ${PARALLELTEST_LIST})
add_custom_target(clean-tests COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Data files
FILE( GLOB DATAFILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/data" "${CMAKE_CURRENT_SOURCE_DIR}/data/*" )

foreach(DATAFILE IN LISTS DATAFILES)
    #add_custom_command( 
    #    TARGET "levelset_001" 
    #    PRE_BUILD
    #    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/${DATAFILE}" "${CMAKE_CURRENT_BINARY_DIR}/data/${DATAFILE}"
    #    ) 
endforeach()

