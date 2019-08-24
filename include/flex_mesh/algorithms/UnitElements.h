#pragma once

#include "flex_mesh/common/Common.h"
#include "flex_mesh/common/Defines.h"
#include "flex_mesh/common/GeometryTypes.h"
#include "flex_mesh/mesh/Mesh.h"
#include "flex_mesh/node/Node.h"
#include "flex_mesh/element/Element.h"

#include <memory>
#include <utility>
#include <vector>

namespace flex_mesh
{
	/// Given a reference to a mesh, adds a linear hex element over the extent of -1 to 1 (3D)
	template<typename MeshType>
	void get_unit_linear_hex(MeshType& mesh);

	/// Returns a mesh (unique_ptr) with a linear hex element over the extent of -1 to 1 (3D)
	template<typename MeshType>
	static std::unique_ptr<MeshType> get_unit_linear_hex_ptr();

	/////////////////
	// Definitions //
	/////////////////

	template<typename MeshType>
	void get_unit_linear_hex(MeshType& mesh)
	{
		auto n1 = mesh.add_node_by_coord({ -1, -1, -1 });
		auto n2 = mesh.add_node_by_coord({ -1, 1, -1 });
		auto n3 = mesh.add_node_by_coord({ -1, 1, 1 });
		auto n4 = mesh.add_node_by_coord({ -1, -1, 1 });
		auto n5 = mesh.add_node_by_coord({ 1, -1, -1 });
		auto n6 = mesh.add_node_by_coord({ 1, 1, -1 });
		auto n7 = mesh.add_node_by_coord({ 1, 1, 1 });
		auto n8 = mesh.add_node_by_coord({ 1, -1, 1 });

		auto e = mesh.elements().add();
		e->initialize(flex_mesh::ElementType::H3D8N, {n1, n2, n3, n4, n5, n6, n7, n8});
		mesh.update_global_numbering();
	}

	template<typename MeshType>
	std::unique_ptr<MeshType> get_unit_linear_hex_ptr()
	{
		std::unique_ptr<MeshType> mesh_ptr = std::make_unique<MeshType>();
		get_unit_linear_hex(*mesh_ptr);
		return std::move(mesh_ptr);
	}
}