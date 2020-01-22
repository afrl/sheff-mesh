cmake_minimum_required(VERSION 3.8)

cmake_policy(SET CMP0048 NEW)
project(FLEXMESH VERSION 1.0.0 LANGUAGES CXX)

# Set package toolchain location for Linux. To enable efficient switching between
# toolchains for cross-compiling, we modify the location where CMake looks for 
# packages using the find_package command by symbolically linking the default location,
# $HOME/.cmake/packages, to a similar location with the folder name given by the 
# compiler family. 
if (UNIX)
    execute_process(COMMAND "mkdir" "-p" "$ENV{HOME}/.cmake/${CMAKE_CXX_COMPILER_ID}"
      RESULT_VARIABLE toolchain_res OUTPUT_VARIABLE toolchain_out)
    if (toolchain_res)
        message(ERROR "Couldn't create toolchain dir: $ENV{HOME}/.cmake/${CMAKE_CXX_COMPILER_ID} ${toolchain_out}")
    endif ()
    execute_process(COMMAND "ln" "-sfn" "$ENV{HOME}/.cmake/${CMAKE_CXX_COMPILER_ID}" "$ENV{HOME}/.cmake/packages"
      RESULT_VARIABLE toolchain_res OUTPUT_VARIABLE toolchain_out)
    if (toolchain_res)
        message(ERROR "Couldn't create toolchain link: $ENV{HOME}/.cmake/packages ${toolchain_out}")
    endif ()
endif()


# Options
option(FLEXMESH_BuildTests "Build the unit tests." OFF)
option(FLEXMESH_USE_FEAIO "Indicates if the build will support the fea_io file format, requiring the library." ON)

# For debugging, print generator, compiler, and options
message(STATUS "Gen: ${CMAKE_GENERATOR}, CXX: ${CMAKE_CXX_COMPILER_ID}, Mode: ${CMAKE_BUILD_TYPE}")
message(STATUS "Tests: ${FLEXMESH_BuildTests}, FEAIO: ${FLEXMESH_USE_FEAIO}")

# Add cmake dir to module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

##################
# Checking input #
##################

include(${CMAKE_SOURCE_DIR}/cmake/CommonInputChecks.cmake)

########################
# Declare dependencies #
########################


find_package(SimpleMath CONFIG REQUIRED HINTS 
  "${CMAKE_SOURCE_DIR}/../simple_math/install_dir" "${CMAKE_SOURCE_DIR}/../../simple_math/install_dir" "${CMAKE_SOURCE_DIR}/modules/simple_math/install_dir")
find_package(Utility CONFIG REQUIRED HINTS 
  "${CMAKE_SOURCE_DIR}/../utility/install_dir" "${CMAKE_SOURCE_DIR}/../../utility/install_dir" "${CMAKE_SOURCE_DIR}/modules/utility/install_dir")
if(FLEXMESH_USE_FEAIO)
find_package(FeaIo CONFIG REQUIRED HINTS 
  "${CMAKE_SOURCE_DIR}/../fea_io/install_dir" "${CMAKE_SOURCE_DIR}/../../fea_io/install_dir" "${CMAKE_SOURCE_DIR}/modules/fea_io/install_dir")
endif()

#######################################
# Setting up FLEXMESH header-only target #
#######################################

# Set the variables for FLEXMESH
set(FLEXMESH_INCLUDE ${CMAKE_SOURCE_DIR}/include)
file (GLOB_RECURSE FLEXMESH_HEADER_LIST "${FLEXMESH_INCLUDE}/*.h")
	
# Create header-only target
add_library(flexmesh INTERFACE)
# Add the header includes and headers to the library
target_include_directories(flexmesh INTERFACE
	$<BUILD_INTERFACE:${FLEXMESH_INCLUDE}>
	$<INSTALL_INTERFACE:include>)
# Sources are only header files, and this is at the interface level so
# downstream targets will have these headers in the VS project
target_sources(flexmesh INTERFACE ${FLEXMESH_HEADER_LIST})
# Add an alias so that a library can be used inside the build tree, eg when testing
add_library(FlexMesh::FlexMesh ALIAS flexmesh)
# Add C++ 14 flag (at interface level since this is a header-only
# library and downstream executables will also need the flag)
target_compile_features(flexmesh INTERFACE cxx_std_14)
# Link dependencies
target_link_libraries(flexmesh INTERFACE SimpleMath::SimpleMath)

if(FLEXMESH_USE_FEAIO AND FeaIo_FOUND)
	target_compile_definitions(flexmesh INTERFACE USE_FEAIO)
	target_link_libraries(flexmesh INTERFACE FeaIo::FeaIo)
endif()

# Assign source groups
include(${CMAKE_SOURCE_DIR}/cmake/Utilities.cmake)
assign_source_group("flex_mesh" ${FLEXMESH_HEADER_LIST})
assign_source_group_SimpleMath()
assign_source_group_Utility()
if(FLEXMESH_USE_FEAIO)
    assign_source_group_FeaIo()
endif()

#################################
# Add unit test targets #
#################################
    
# Build tests
if(FLEXMESH_BuildTests)
    include(CTest)
    enable_testing()
    add_subdirectory(unit_tests)
endif()

################
# Installation #
################

include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/FlexMesh)

# Install the header files
install(DIRECTORY ${FLEXMESH_INCLUDE}/
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
# This is required so that the exported target has the name FlexMesh and not flexmesh
set_target_properties(flexmesh PROPERTIES EXPORT_NAME FlexMesh)

# Create a ConfigVersion.cmake file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
	${CMAKE_CURRENT_BINARY_DIR}/FlexMeshConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

# Create the package config file for find_package(FlexMesh)
configure_package_config_file(
	${CMAKE_CURRENT_LIST_DIR}/cmake/FlexMeshConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/FlexMeshConfig.cmake
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

# Install the config and config version files
# Note: if you write custom find modules, add them here
# Refer to https://github.com/pabloariasal/modern-cmake-sample/blob/cfda7b51ff915409fd61d4f2e75b63032bedb182/libjsonutils/CMakeLists.txt
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/FlexMeshConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/FlexMeshConfigVersion.cmake
	${CMAKE_CURRENT_BINARY_DIR}/FlexMeshTargets.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)

###########
# Exports #
###########

# Export just the flexmesh target (not tests, etc.)
export(
    TARGETS flexmesh
	FILE ${CMAKE_CURRENT_BINARY_DIR}/FlexMeshTargets.cmake
    NAMESPACE FlexMesh::
)

# Register package in user's package registry
export(PACKAGE FlexMesh)
