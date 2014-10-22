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

#Build test_3D_edo.cpp
SET(test_3D_edo_src test_3D_edo.cpp)

add_executable(test_3D_edo ${test_3D_edo_src})

target_link_libraries(test_3D_edo mpi)
TARGET_LINK_LIBRARIES(test_3D_edo pablo)

#Build test_2D3D_debug.cpp
SET(test_2D3D_debug_src test_2D3D_debug.cpp)

add_executable(test_2D3D_debug ${test_2D3D_debug_src})

target_link_libraries(test_2D3D_debug mpi)
TARGET_LINK_LIBRARIES(test_2D3D_debug pablo)

#Build test_2D_restart_debug.cpp
SET(test_2D_restart_debug_src test_2D_restart_debug.cpp)
add_executable(test_2D_restart_debug ${test_2D_restart_debug_src})

target_link_libraries(test_2D_restart_debug mpi)
TARGET_LINK_LIBRARIES(test_2D_restart_debug pablo)
