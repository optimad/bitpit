#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

#Name your project here
project(PABLO)

ADD_DEFINITIONS(-std=c++11)
SET(WITHOUT_MPI 0 CACHE BOOL "Set WITHOUT_MPI to 1 if you don't have an installed implementation of MPI or you don't want to use a parallel PABLO")
SET(COMPILER "gcc" CACHE STRING "Set COMIPLER to 'intel' if you want to use the Intel Compiler")
IF(WITHOUT_MPI EQUAL 1)
message("-- PABLO Serial Version")
ADD_DEFINITIONS(-DNOMPI=1)
ENDIF(WITHOUT_MPI EQUAL 1)
IF(WITHOUT_MPI EQUAL 0)
message("-- PABLO Parallel Version")
ENDIF(WITHOUT_MPI EQUAL 0)

include_directories(include)

IF(WITHOUT_MPI EQUAL 0)

IF(COMPILER EQUAL intel)
SET(CMAKE_C_COMPILER mpiicc)
SET(CMAKE_CXX_COMPILER mpicpc)
ELSE(COMPILER EQUAL "intel")
SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpic++)
ENDIF(COMPILER EQUAL "intel")

ELSE(WITHOUT_MPI EQUAL 0)

IF(COMPILER EQUAL intel)
SET(CMAKE_C_COMPILER icc)
SET(CMAKE_CXX_COMPILER icpc)
ELSE(COMPILER EQUAL "intel")
SET(CMAKE_C_COMPILER gcc)
SET(CMAKE_CXX_COMPILER g++)
ENDIF(COMPILER EQUAL intel)

ENDIF(WITHOUT_MPI EQUAL 0)

MESSAGE("Compiling with ${COMPILER}")

SET(DEBUG 0 CACHE STRING "Set DEBUG = 1 if you want debug flags to be passed to the compiler. Otherwise, PABLO will be optimized")

IF(${DEBUG})
IF(COMPILER EQUAL "intel")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")
ELSE(COMPILER EQUAL "intel")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ggdb -O0 -fmessage-length=0 -fno-weak")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb -O0 -fmessage-length=0 -fno-weak")
ENDIF(COMPILER EQUAL "intel")
ELSE(${DEBUG})
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif(${DEBUG})

ADD_SUBDIRECTORY(lib)

SET(COMPILE_TESTS 1 CACHE STRING "Set COMPILE_TESTS = 0 if you dont't want to compile PABLO's tests. Otherwise, PABLO's tests will be compiled")
IF(${COMPILE_TESTS})
ADD_SUBDIRECTORY(test)
ENDIF(${COMPILE_TESTS})

SET(CMAKE_INSTALL_PREFIX:PATH /usr/local/)
INSTALL(DIRECTORY include/ DESTINATION ${PROJECT_NAME}/include)
EXEC_PROGRAM(pwd . OUTPUT_VARIABLE BUILD_DIR)
INSTALL(FILES ${BUILD_DIR}/lib/libPABLO.a DESTINATION lib)

UNSET(DEBUG CACHE)
UNSET(COMPILE_TESTS CACHE)
