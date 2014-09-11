#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

#Name your project here
project(OPTIMAD_base_library)

ADD_DEFINITIONS(-std=c++0x)

include_directories(include)

ADD_SUBDIRECTORY(lib )
ADD_SUBDIRECTORY(test_lib)
