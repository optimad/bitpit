file(GLOB PABLO_src "*.cpp")
#SET(PABLO_src {*.cpp})
add_library(PABLO ${PABLO_src})
target_link_libraries(PABLO mpi)