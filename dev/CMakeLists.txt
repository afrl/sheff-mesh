set(FLEX_MESH_DEV_TARGET_NAME ${PROJECT_NAME}_dev)
add_executable(${FLEX_MESH_DEV_TARGET_NAME} source/main.cpp)
target_link_libraries(${FLEX_MESH_DEV_TARGET_NAME} PRIVATE ${FLEX_MESH_TARGET_NAME})
assign_source_group(${FLEX_MESH_HEADER_LIST})
