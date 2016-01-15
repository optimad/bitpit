# Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

#------------------------------------------------------------------------------------#
# Subdirectories
#------------------------------------------------------------------------------------#

if (DEFINED TEST_LIST)
	unset(TEST_LIST)
endif ()

# Modules
foreach(MODULE_NAME IN LISTS BITPIT_MODULE_LIST)
	isModuledEnabled(${MODULE_NAME} MODULE_ENABLED)
	if (MODULE_ENABLED)
		if (IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE_NAME}")
			string(TOUPPER ${MODULE_NAME} UPPERCASE_MODULE_NAME)
			add_subdirectory(${MODULE_NAME})
			list(APPEND TEST_LIST "${${UPPERCASE_MODULE_NAME}_TEST_LIST}")
		endif ()
	endif ()
endforeach()

#------------------------------------------------------------------------------------#
# Targets
#------------------------------------------------------------------------------------#

add_custom_target(tests DEPENDS ${TEST_LIST})
add_custom_target(clean-tests COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
