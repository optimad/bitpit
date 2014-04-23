#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

#Name your project here

#Build test_pablo.cpp
#SET(test_pablo_src test_pablo.cpp)
#
#add_executable(test_pablo ${test_pablo_src})
#
#target_link_libraries(test_pablo mpi)
#TARGET_LINK_LIBRARIES(test_pablo pablo)
#
#Build test_marco.cpp
#SET(test_marco_src test_marco.cpp)
#
#add_executable(test_marco ${test_marco_src})
#
#target_link_libraries(test_marco mpi)
#TARGET_LINK_LIBRARIES(test_marco pablo)
#
#Build test_edo.cpp
#SET(test_edo_src test_edo.cpp)
#
#add_executable(test_edo ${test_edo_src})
#
#target_link_libraries(test_edo mpi)
#TARGET_LINK_LIBRARIES(test_edo pablo)

#Build test_2D3D.cpp
SET(test_2D3D_src test_2D3D.cpp)

add_executable(test_2D3D ${test_2D3D_src})

target_link_libraries(test_2D3D mpi)
TARGET_LINK_LIBRARIES(test_2D3D pablo)

#Build test_2D3D_edo.cpp
SET(test_2D3D_src test_2D3D_edo.cpp)

add_executable(test_2D3D_edo ${test_2D3D_src})

target_link_libraries(test_2D3D_edo mpi)
TARGET_LINK_LIBRARIES(test_2D3D_edo pablo)
