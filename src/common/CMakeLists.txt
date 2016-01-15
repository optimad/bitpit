# CMake settings
cmake_minimum_required(VERSION 2.8)

# Reset the variable
set(HEADERS "")

# Generate the reference header files
set(HEADER_TEMPLATE "BitPit_template.hpp.in")
foreach (GROUP_NAME IN LISTS BITPIT_GROUP_LIST)
	set(REFERENCE_HEADER "BitPit_${GROUP_NAME}.hpp")
	string(TOUPPER ${GROUP_NAME} UPPER_GROUP_NAME)

	set(INCLUDE_LIST "")
	foreach(MODULE_NAME IN LISTS BITPIT_${UPPER_GROUP_NAME}_MODULE_LIST)
		isModuledEnabled(${MODULE_NAME} MODULE_ENABLED)
		if (MODULE_ENABLED)
			set(INCLUDE_LIST "${INCLUDE_LIST}#include <BitPit_Base_${MODULE_NAME}.hpp>\n")
		endif ()
	endforeach ()

	CONFIGURE_FILE("${HEADER_TEMPLATE}" "${REFERENCE_HEADER}")
	list(APPEND HEADERS "${CMAKE_CURRENT_BINARY_DIR}/${REFERENCE_HEADER}")

	unset(INCLUDE_LIST)
	unset(REFERENCE_HEADER)
endforeach ()

# Add common files to targets
file(GLOB HEADER_FILES "*.hpp" "*.tpp")
list(APPEND HEADERS "${HEADER_FILES}")
unset(HEADER_FILES)

set(COMMON_HEADERS "${HEADERS}" CACHE INTERNAL "Common headers" FORCE)
install(FILES ${COMMON_HEADERS} DESTINATION include/${CMAKE_PROJECT_NAME})
