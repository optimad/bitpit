# Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

string(TOUPPER ${CMAKE_PROJECT_NAME} FIND_MODULE_FILENAME)
set(FIND_MODULE_FILENAME "Find${FIND_MODULE_FILENAME}.cmake")

add_custom_target(clean-cmake COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

install(FILES "${FIND_MODULE_FILENAME}" DESTINATION lib/${CMAKE_PROJECT_NAME})
install(FILES "LibFindMacros.cmake" DESTINATION lib/${CMAKE_PROJECT_NAME})