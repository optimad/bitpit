#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

# List of tests
set(SERIAL_TESTS "")

list(APPEND TESTS "pablo_001")
list(APPEND TESTS "pablo_002")
list(APPEND TESTS "pablo_003")
list(APPEND TESTS "pablo_004")
if (NOT ONLY_PABLO)
    list(APPEND TESTS "ucartmesh_001")
    list(APPEND TESTS "ucartmesh_002")
    list(APPEND TESTS "surftri_001")
    list(APPEND TESTS "surftri_002")
    list(APPEND TESTS "surftri_003")
    list(APPEND TESTS "surftri_004")
    list(APPEND TESTS "surftri_005")
    list(APPEND TESTS "voltri_001")
    list(APPEND TESTS "patchman_001")
    list(APPEND TESTS "patchman_002")
    list(APPEND TESTS "patchman_003")
endif()

if (ENABLE_MPI)
    set(PARALLEL_TEST "")
    list(APPEND PARALLEL_TESTS "parallel_pablo_001")
endif()

set(TEST_LIST "${TESTS}" CACHE INTERNAL "List of serial tests" FORCE)
set(PARALLELTEST_LIST "${PARALLEL_TESTS}" CACHE INTERNAL "List of parallel tests" FORCE)


include_directories("${PROJECT_SOURCE_DIR}/src/common")
include_directories("${PROJECT_SOURCE_DIR}/src/pablo")
if (NOT ONLY_PABLO)
    include_directories("${PROJECT_SOURCE_DIR}/src/ucartmesh")
    include_directories("${PROJECT_SOURCE_DIR}/src/surftri")
    include_directories("${PROJECT_SOURCE_DIR}/src/voltri")
    include_directories("${PROJECT_SOURCE_DIR}/src/patchman")
endif()

include_directories(${BITPIT_INCLUDE_DIRS})

#Rules to build the tests
foreach(TEST_NAME IN LISTS TEST_LIST PARALLELTEST_LIST)
    set(TEST_SOURCES "${TEST_NAME}.cpp")

    add_executable(${TEST_NAME} "${TEST_SOURCES}")
    target_link_libraries(${TEST_NAME} ${BITP_MESH_LIBRARY})
    target_link_libraries(${TEST_NAME} ${BITPIT_LIBRARY})
endforeach()

add_custom_target(tests DEPENDS ${TEST_LIST} ${PARALLELTEST_LIST})
add_custom_target(clean-tests COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Data files
FILE( GLOB DATAFILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/data" "${CMAKE_CURRENT_SOURCE_DIR}/data/*" )

foreach(DATAFILE IN LISTS DATAFILES)
if (NOT ONLY_PABLO)
    add_custom_command( 
        TARGET "surftri_001" 
        PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/${DATAFILE}" "${CMAKE_CURRENT_BINARY_DIR}/data/${DATAFILE}"
        ) 
endif()
endforeach()

