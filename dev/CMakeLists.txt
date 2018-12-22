add_executable(flex_mesh_dev source/main.cpp)
target_link_libraries(flex_mesh_dev PRIVATE FlexMesh::FlexMesh)

include(${CMAKE_SOURCE_DIR}/cmake/Utilities.cmake)
assign_source_group("flex_mesh" ${FLEXMESH_HEADER_LIST})
if(FLEXMESH_USE_FEAIO)
	assign_source_group_feaio()
endif()
