#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

#Name your project here

#Build test0.cpp
SET(test0_src test0.cpp)

add_executable(test0 ${test0_src})

target_link_libraries(test0 mpi)
TARGET_LINK_LIBRARIES(test0 pablo)

#Build test1.cpp
SET(test1_src test1.cpp)

add_executable(test1 ${test1_src})

target_link_libraries(test1 mpi)
TARGET_LINK_LIBRARIES(test1 pablo)

#Build test2.cpp
SET(test2_src test2.cpp)

add_executable(test2 ${test2_src})

target_link_libraries(test2 mpi)
TARGET_LINK_LIBRARIES(test2 pablo)

#Build test3.cpp
SET(test3_src test3.cpp)

add_executable(test3 ${test3_src})

target_link_libraries(test3 mpi)
TARGET_LINK_LIBRARIES(test3 pablo)
