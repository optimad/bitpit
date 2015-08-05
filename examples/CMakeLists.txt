#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

# Add a target to generate the examples
option(BUILD_EXAMPLES "Create and install the examples" OFF)

if(BUILD_EXAMPLES)
	# List of examples
	set(EXAMPLE_LIST "")
	list(APPEND EXAMPLE_LIST "createCartesianMesh")

	#Rules to build the examples
	foreach(EXAMPLE_NAME IN LISTS EXAMPLE_LIST)
		set(EXAMPLE_SOURCES "${EXAMPLE_NAME}.cpp")

		add_executable(${EXAMPLE_NAME} "${EXAMPLE_SOURCES}")
		target_link_libraries(${EXAMPLE_NAME} ${PATCHMAN_LIBRARY})
		if (PABLO_FOUND)
			target_link_libraries(${EXAMPLE_NAME} ${PABLO_LIBRARIES})
		endif()
		target_link_libraries(${EXAMPLE_NAME} ${VTK_LIBRARIES})
		install(TARGETS ${EXAMPLE_NAME} DESTINATION doc/${CMAKE_PROJECT_NAME}-${PATCHMAN_VERSION}/examples)
	endforeach()

	add_custom_target(examples DEPENDS ${EXAMPLE_LIST})
endif()