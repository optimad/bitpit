file(GLOB pablo_src "*.cpp")
#SET(pablo_src {*.cpp})
add_library(pablo ${pablo_src})
target_link_libraries(pablo mpi)