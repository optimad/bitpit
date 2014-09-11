#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

#Name your project here

file(GLOB test_obl_src "*.cpp")
add_executable(test_obl ${test_obl_src})


#TARGET_LINK_LIBRARIES(test_obl obl)
