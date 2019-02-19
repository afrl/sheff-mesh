#pragma once

#include "flex_mesh/common/Common.h"
#include "flex_mesh/common/Defines.h"
#include "flex_mesh/common/GeometryTypes.h"
#include "flex_mesh/mesh/Mesh.h"
#include "flex_mesh/node/Node.h"
#include "flex_mesh/element/Element.h"

#include <utility>
#include <vector>

namespace flex_mesh
{
	/// Returns the nodes and elements respectively for the given partition in the mesh
	template<typename MeshType>
	static MeshType get_unit_linear_hex();

	/////////////////
	// Definitions //
	/////////////////

	template<typename MeshType>
	MeshType get_unit_linear_hex()
	{
		MeshType mesh;
		auto n1 = mesh.add_node_by_coord({ -1, -1, -1 });
		auto n2 = mesh.add_node_by_coord({ 1, -1, -1 });
		auto n3 = mesh.add_node_by_coord({ 1, 1, -1 });
		auto n4 = mesh.add_node_by_coord({ -1, 1, -1 });
		auto n5 = mesh.add_node_by_coord({ -1, -1, 1 });
		auto n6 = mesh.add_node_by_coord({ 1, -1, 1 });
		auto n7 = mesh.add_node_by_coord({ 1, 1, 1 });
		auto n8 = mesh.add_node_by_coord({ -1, 1, 1 });

		auto e = mesh.elements().add();
		e->initialize(flex_mesh::ElementType::H3D8N, {n1, n2, n3, n4, n5, n6, n7, n8});
		return mesh;
	}
}