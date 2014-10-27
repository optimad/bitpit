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

#Build test4.cpp
SET(test4_src test4.cpp)

add_executable(test4 ${test4_src})

target_link_libraries(test4 mpi)
TARGET_LINK_LIBRARIES(test4 pablo)

#Build test5.cpp
SET(test5_src test5.cpp)

add_executable(test5 ${test5_src})

target_link_libraries(test5 mpi)
TARGET_LINK_LIBRARIES(test5 pablo)


#Build test12.cpp
SET(test12_src test12.cpp)

add_executable(test12 ${test12_src})

target_link_libraries(test12 mpi)
TARGET_LINK_LIBRARIES(test12 pablo)

#Build test13.cpp
SET(test13_src test13.cpp)

add_executable(test13 ${test13_src})

target_link_libraries(test13 mpi)
TARGET_LINK_LIBRARIES(test13 pablo)

#Build test14.cpp
SET(test14_src test14.cpp)

add_executable(test14 ${test14_src})

target_link_libraries(test14 mpi)
TARGET_LINK_LIBRARIES(test14 pablo)

#Build test15.cpp
SET(test15_src test15.cpp)

add_executable(test15 ${test15_src})

target_link_libraries(test15 mpi)
TARGET_LINK_LIBRARIES(test15 pablo)

#Build test16.cpp
SET(test16_src test16.cpp)

add_executable(test16 ${test16_src})

target_link_libraries(test16 mpi)
TARGET_LINK_LIBRARIES(test16 pablo)

#Build test104.cpp
SET(test104_src test104.cpp)

add_executable(test104 ${test104_src})

target_link_libraries(test104 mpi)
TARGET_LINK_LIBRARIES(test104 pablo)

#Build test116.cpp
SET(test116_src test116.cpp)

add_executable(test116 ${test116_src})

target_link_libraries(test116 mpi)
TARGET_LINK_LIBRARIES(test116 pablo)
