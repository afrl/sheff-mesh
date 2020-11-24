cmake_minimum_required(VERSION 3.8)

cmake_policy(SET CMP0048 NEW)
project(MESH VERSION 1.0.0 LANGUAGES CXX)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # require at least gcc 6.0
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.0)
        message( FATAL_ERROR "GCC version must be at least 6.0! (Found: ${CMAKE_CXX_COMPILER_VERSION})" )
    endif()
endif()

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
option(MESH_BuildTests "Build the unit tests." OFF)
option(MESH_USE_FEAIO "Indicates if the build will support the fea_io file format, requiring the library." ON)

# For debugging, print generator, compiler, and options
message(STATUS "Gen: ${CMAKE_GENERATOR}, CXX: ${CMAKE_CXX_COMPILER_ID}, Mode: ${CMAKE_BUILD_TYPE}")
message(STATUS "Tests: ${MESH_BuildTests}, FEAIO: ${MESH_USE_FEAIO}")

# Add cmake dir to module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

##################
# Checking input #
##################

include(${CMAKE_SOURCE_DIR}/cmake/CommonInputChecks.cmake)

########################
# Declare dependencies #
########################

find_package(RapidJSON 1.1.0 CONFIG REQUIRED HINTS 
  "${CMAKE_SOURCE_DIR}/../third_party/rapidjson_install_dir" "${CMAKE_SOURCE_DIR}/../../third_party/rapidjson_install_dir" "${CMAKE_SOURCE_DIR}/modules/third_party/rapidjson_install_dir")
find_package(SecondPartyMetis CONFIG REQUIRED HINTS 
  "${CMAKE_SOURCE_DIR}/../second_party/second_party_metis/install_dir" "${CMAKE_SOURCE_DIR}/../second_party_metis/install_dir" "${CMAKE_SOURCE_DIR}/modules/second_party/second_party_metis/install_dir")
  find_package(SecondPartyMesh CONFIG REQUIRED HINTS 
  "${CMAKE_SOURCE_DIR}/../second_party/second_party_mesh/install_dir" "${CMAKE_SOURCE_DIR}/../second_party_mesh/install_dir" "${CMAKE_SOURCE_DIR}/modules/second_party/second_party_mesh/install_dir")
find_package(SimpleMath CONFIG REQUIRED HINTS 
  "${CMAKE_SOURCE_DIR}/../simple_math/install_dir" "${CMAKE_SOURCE_DIR}/../../simple_math/install_dir" "${CMAKE_SOURCE_DIR}/modules/simple_math/install_dir")
find_package(Utility CONFIG REQUIRED HINTS 
  "${CMAKE_SOURCE_DIR}/../utility/install_dir" "${CMAKE_SOURCE_DIR}/../../utility/install_dir" "${CMAKE_SOURCE_DIR}/modules/utility/install_dir")

if(MESH_USE_FEAIO)
  find_package(FeaIo CONFIG REQUIRED HINTS 
    "${CMAKE_SOURCE_DIR}/../fea_io/install_dir" "${CMAKE_SOURCE_DIR}/../../fea_io/install_dir" "${CMAKE_SOURCE_DIR}/modules/fea_io/install_dir")
endif()

######################################
# Setting up MESH_INCLUDE target #
######################################

# Set the variables for MESH
set( MESH_INCLUDE     ${CMAKE_SOURCE_DIR}/include )
set( MESH_SRC         ${CMAKE_SOURCE_DIR}/source )

file( GLOB_RECURSE MESH_HEADER_LIST "${MESH_INCLUDE}/*.h" )
file( GLOB_RECURSE MESH_SOURCE_LIST "${MESH_SRC}/*.h" "${MESH_SRC}/*.cpp" )
	
# Create header-only target
add_library(mesh STATIC)
# Add the header includes and headers to the library
target_include_directories(mesh PUBLIC
    $<BUILD_INTERFACE:${MESH_INCLUDE}>
    $<INSTALL_INTERFACE:include>)
# Add headers in the include folder as PUBLIC so downstream
# targets have access to them
target_sources(mesh PUBLIC ${MESH_HEADER_LIST})
# Add headers and sources in the source folder as PRIVATE as they
# will not be needed downstream
target_sources(mesh PRIVATE ${MESH_SOURCE_LIST})
# Add an alias so that a library can be used inside the build tree, eg when testing
add_library(Mesh::Mesh ALIAS mesh)
# Add C++ 14 flag
target_compile_features(mesh PUBLIC cxx_std_14)

# Link dependencies
target_include_directories(mesh PUBLIC ${RapidJSON_INCLUDE_DIR})
target_link_libraries(mesh PUBLIC SecondPartyMesh::SecondPartyMesh)
target_link_libraries(mesh PUBLIC SecondPartyMetis::SecondPartyMetis)
target_link_libraries(mesh PUBLIC SimpleMath::SimpleMath)

if(MESH_USE_FEAIO AND FeaIo_FOUND)
	target_compile_definitions(mesh PUBLIC USE_FEAIO)
	target_link_libraries(mesh PUBLIC FeaIo::FeaIo)
endif()

# Assign source groups
include(${CMAKE_SOURCE_DIR}/cmake/Utilities.cmake)
assign_source_group("mesh" ${MESH_HEADER_LIST})
assign_source_group("mesh" ${MESH_SOURCE_LIST})
assign_source_group_SecondPartyMesh()
assign_source_group_SecondPartyMetis()
assign_source_group_SimpleMath()
assign_source_group_Utility()
if(MESH_USE_FEAIO)
    assign_source_group_FeaIo()
endif()

#################################
# Add unit test targets #
#################################
    
# Build tests
if(MESH_BuildTests)
    include(CTest)
    enable_testing()
    add_subdirectory(unit_tests)
endif()

################
# Installation #
################

include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/Mesh)

# Install the mesh target (library and headers)
install(
    TARGETS mesh
    EXPORT MeshTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    INCLUDES DESTINATION include
)
# This is required so that the exported target has the name Mesh and not mesh
set_target_properties(mesh PROPERTIES EXPORT_NAME Mesh)

# Create a ConfigVersion.cmake file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/MeshConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

# Create the package config file for find_package(Mesh)
configure_package_config_file(
    ${CMAKE_CURRENT_LIST_DIR}/cmake/MeshConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/MeshConfig.cmake
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

# Install the config and config version files
# Note: if you write custom find modules, add them here
# Refer to https://github.com/pabloariasal/modern-cmake-sample/blob/cfda7b51ff915409fd61d4f2e75b63032bedb182/libjsonutils/CMakeLists.txt
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/MeshConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/MeshConfigVersion.cmake
	  ${CMAKE_CURRENT_BINARY_DIR}/MeshTargets.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)

###########
# Exports #
###########

# Export just the mesh target (not tests, etc.)
export(
    TARGETS mesh
    FILE ${CMAKE_CURRENT_BINARY_DIR}/MeshTargets.cmake
    NAMESPACE Mesh::
)

# Register package in user's package registry
export(PACKAGE Mesh)
