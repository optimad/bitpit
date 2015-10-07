#------------------------------------------------------------------------------------#
# Subdirectories
#------------------------------------------------------------------------------------#

# Common
add_subdirectory(common)

# Base
add_subdirectory(base)

# IO
add_subdirectory(IO)

# SA
add_subdirectory(SA)

# Library
add_custom_target(library DEPENDS ${BASE_LIBRARY})
add_custom_target(clean-library COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
