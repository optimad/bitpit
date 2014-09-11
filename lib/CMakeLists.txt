#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

#Name your project here

file(GLOB lib_base_src "*.cpp")

add_library(opt_base STATIC ${lib_base_src})
