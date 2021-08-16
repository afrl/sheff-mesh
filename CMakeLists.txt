cmake_minimum_required(VERSION 3.8)

cmake_policy(SET CMP0048 NEW)
project(MESH VERSION 1.0.0 LANGUAGES CXX)

# Options
option(SHEFF_BuildTests "Build the unit tests." OFF)
option(SHEFF_USE_FEAIO "Indicates if the build will support the fea_io file format, requiring the library." ON)

# For debugging, print generator, compiler, and options
message(STATUS "Gen: ${CMAKE_GENERATOR}, CXX: ${CMAKE_CXX_COMPILER_ID}, Mode: ${CMAKE_BUILD_TYPE}")
message(STATUS "Tests: ${SHEFF_BuildTests}, FEAIO: ${SHEFF_USE_FEAIO}")

# Add cmake dir to module path
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

##################
# Checking input #
##################

include(${PROJECT_SOURCE_DIR}/cmake/CommonInputChecks.cmake)

########################
# Declare dependencies #
########################

if(NOT RapidJSON_FOUND)
    message(STATUS "SecondPartyMath did not already exist, finding package.")
    find_package(RapidJSON 1.1.0 CONFIG REQUIRED HINTS
                 "${PROJECT_SOURCE_DIR}/../../third_party/rapidjson_install_dir")
endif()
if(NOT TARGET Utility::Utility)
    message(STATUS "Utility did not already exist, finding package.")
    find_package(Utility CONFIG REQUIRED HINTS "${PROJECT_SOURCE_DIR}/../../first_party/utility/build")
endif()
if(NOT TARGET SecondPartyMesh::SecondPartyMesh)
    message(STATUS "SecondPartyMesh did not already exist, finding package.")
    find_package(SecondPartyMesh CONFIG REQUIRED HINTS "${PROJECT_SOURCE_DIR}/../../second_party/second_party_mesh/build")
endif()
if(NOT TARGET SecondPartyMetis::SecondPartyMetis)
    message(STATUS "SecondPartyMetis did not already exist, finding package.")
    find_package(SecondPartyMetis CONFIG REQUIRED HINTS "${PROJECT_SOURCE_DIR}/../../second_party/second_party_metis/build")
endif()
if(NOT TARGET SimpleMath::SimpleMath)
    message(STATUS "SimpleMath did not already exist, finding package.")
    find_package(SimpleMath CONFIG REQUIRED HINTS "${PROJECT_SOURCE_DIR}/../../first_party/simple_math/build")
endif()

if(SHEFF_USE_FEAIO)
    if(NOT TARGET FeaIo::FeaIo)
        message(STATUS "FeaIo did not already exist, finding package.")
        find_package(FeaIo CONFIG REQUIRED HINTS "${PROJECT_SOURCE_DIR}/../../first_party/fea_io/build")
    endif()
endif()

######################################
# Setting up MESH_INCLUDE target #
######################################

# Set the variables for MESH
set( MESH_INCLUDE     ${PROJECT_SOURCE_DIR}/include )
set( MESH_SRC         ${PROJECT_SOURCE_DIR}/source )

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
set_target_properties(mesh PROPERTIES FOLDER "sheff")
# Link dependencies
target_include_directories(mesh PUBLIC ${RapidJSON_INCLUDE_DIR})
target_link_libraries(mesh PUBLIC SecondPartyMesh::SecondPartyMesh)
target_link_libraries(mesh PUBLIC SecondPartyMetis::SecondPartyMetis)
target_link_libraries(mesh PUBLIC SimpleMath::SimpleMath)

if(SHEFF_USE_FEAIO)
	target_compile_definitions(mesh PUBLIC USE_FEAIO)
	target_link_libraries(mesh PUBLIC FeaIo::FeaIo)
endif()

###########################
# Assigning source groups #
###########################

include("$ENV{SHEFF_TOOLS}/cmake/IncludeSHEFF.cmake")
register_sources_sheff("mesh"
                       "${MESH_HEADER_LIST}"
                       "${MESH_SOURCE_LIST}")
assign_source_groups_sheff()

#################################
# Add unit test targets #
#################################
    
# Build tests
if(SHEFF_BuildTests)
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
