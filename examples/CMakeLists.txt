#---------------------------------------------------------------------------
#
#  bitpit
#
#  Copyright (C) 2015-2021 OPTIMAD engineering Srl
#
#  -------------------------------------------------------------------------
#  License
#  This file is part of bitpit.
#
#  bitpit is free software: you can redistribute it and/or modify it
#  under the terms of the GNU Lesser General Public License v3 (LGPL)
#  as published by the Free Software Foundation.
#
#  bitpit is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
#
#---------------------------------------------------------------------------*/

# Examples can be compiled only if all modules are enabled
set(MODULE_LIST "${BITPIT_MODULE_LIST}")
if (NOT BITPIT_ENABLE_MPI)
    list(REMOVE_ITEM MODULE_LIST "${MPI_ONLY_MODULE_LIST}")
endif ()

areModulesEnabled("${MODULE_LIST}" ALL_MODULES_ENABLED)
if (ALL_MODULES_ENABLED)
    option(BITPIT_BUILD_EXAMPLES "Create and install the examples" OFF)
else ()
    unset(BITPIT_BUILD_EXAMPLES CACHE)
    return()
endif ()

if(BITPIT_BUILD_EXAMPLES)
    # List of examples
    set(EXAMPLE_LIST "")
    list(APPEND EXAMPLE_LIST "LA_example_00001")
    list(APPEND EXAMPLE_LIST "LA_example_00002")
    list(APPEND EXAMPLE_LIST "operators_example")
    list(APPEND EXAMPLE_LIST "PABLO_bubbles_2D")
    list(APPEND EXAMPLE_LIST "PABLO_bubbles_3D")
    list(APPEND EXAMPLE_LIST "PABLO_example_00001")
    list(APPEND EXAMPLE_LIST "PABLO_example_00002")
    list(APPEND EXAMPLE_LIST "PABLO_example_00003")
    list(APPEND EXAMPLE_LIST "PABLO_example_00004")
    list(APPEND EXAMPLE_LIST "PABLO_example_00005")
    list(APPEND EXAMPLE_LIST "PABLO_example_00006")
    list(APPEND EXAMPLE_LIST "PABLO_example_00007")
    list(APPEND EXAMPLE_LIST "PABLO_example_00008")
    list(APPEND EXAMPLE_LIST "PABLO_example_00009")
    list(APPEND EXAMPLE_LIST "PABLO_example_00010")
    list(APPEND EXAMPLE_LIST "PABLO_example_00011")
    list(APPEND EXAMPLE_LIST "patchkernel_example_00001")
    list(APPEND EXAMPLE_LIST "volcartesian_example_00001")
    list(APPEND EXAMPLE_LIST "POD_example_00001")
    list(APPEND EXAMPLE_LIST "POD_example_00002")
    list(APPEND EXAMPLE_LIST "POD_example_00003")
    list(APPEND EXAMPLE_LIST "POD_example_00004")
    list(APPEND EXAMPLE_LIST "POD_example_00005")
    list(APPEND EXAMPLE_LIST "POD_application_example_00001")
    list(APPEND EXAMPLE_LIST "RBF_example_00001")
    list(APPEND EXAMPLE_LIST "voloctree_adaptation_example_00001")
    list(APPEND EXAMPLE_LIST "voloctree_mapper_example_00001")
    list(APPEND EXAMPLE_LIST "voloctree_mapper_example_00002")
    list(APPEND EXAMPLE_LIST "voloctree_mapper_example_00003")
    list(APPEND EXAMPLE_LIST "voloctree_mapper_example_00004")

    #Rules to build the examples
    foreach(EXAMPLE_NAME IN LISTS EXAMPLE_LIST)
        set(EXAMPLE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/${EXAMPLE_NAME}.cpp")
        list(APPEND EXAMPLES_SOURCES ${EXAMPLE_SOURCES})

        add_executable(${EXAMPLE_NAME} "${EXAMPLE_SOURCES}")

        target_compile_features(${EXAMPLE_NAME} PUBLIC cxx_std_11)
        set_target_properties(${EXAMPLE_NAME} PROPERTIES CXX_STANDARD 11)
        set_target_properties(${EXAMPLE_NAME} PROPERTIES CXX_STANDARD_REQUIRED ON)

        target_link_libraries(${EXAMPLE_NAME} ${BITPIT_LIBRARY})
        target_link_libraries(${EXAMPLE_NAME} ${BITPIT_EXTERNAL_LIBRARIES})

        if (MSVC)
            target_compile_definitions(${EXAMPLE_NAME} PRIVATE "_CRT_SECURE_NO_WARNINGS")
        endif()

        install(TARGETS ${EXAMPLE_NAME} DESTINATION doc/${PROJECT_NAME}-${BITPIT_VERSION}/examples)

        foreach (MODULE_NAME IN LISTS BITPIT_MODULE_LIST)
            isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)
            if (MODULE_ENABLED)
                addModuleIncludeDirectories(${EXAMPLE_NAME} ${MODULE_NAME})
            endif()
        endforeach ()
    endforeach()

    add_custom_target(examples DEPENDS ${EXAMPLE_LIST})
    add_custom_target(clean-examples COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

    # Data files
    add_custom_command(
        TARGET "POD_application_example_00001" "POD_example_00001" "POD_example_00002" "POD_example_00003" PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/data" "${CMAKE_CURRENT_BINARY_DIR}/data"
        COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_BINARY_DIR}/pod"
    )

    add_custom_command(
        TARGET "RBF_example_00001" PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/RBF_example_00001.param" "${CMAKE_CURRENT_BINARY_DIR}/RBF_example_00001.param"
        COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/RBF_example_00001.stl" "${CMAKE_CURRENT_BINARY_DIR}/RBF_example_00001.stl"
    )

    install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/data" DESTINATION doc/${PROJECT_NAME}-${BITPIT_VERSION}/examples)

    # List of example sources
    set (EXAMPLES_SOURCES ${EXAMPLES_SOURCES} PARENT_SCOPE)
endif()
