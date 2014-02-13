#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

#Name your project here
project(PABLO)

ADD_DEFINITIONS(-std=c++11)

include_directories(include)
SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpic++)

ADD_SUBDIRECTORY(lib)
ADD_SUBDIRECTORY(test)

IF(CMAKE_COMPILER_IS_GNUCC)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0")
ENDIF(CMAKE_COMPILER_IS_GNUCC)
IF(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)


