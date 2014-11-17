#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

#Name your project here
project(PABLO)

ADD_DEFINITIONS(-std=c++11)

include_directories(include)
SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpic++)
#SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ggdb -O0 -fmessage-length=0")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb -O0 -fmessage-length=0")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

ADD_SUBDIRECTORY(lib)
ADD_SUBDIRECTORY(test)


