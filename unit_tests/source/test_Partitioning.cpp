#include <chrono>
#include <functional>
#include <iostream>
#include <iomanip>
#include <vector>

#include "flex_mesh/algorithms/Creation.h"
#include "flex_mesh/algorithms/IO.h"
#include "flex_mesh/algorithms/Partitioning.h"
#include "flex_mesh/node/Node.h"
#include "flex_mesh/element/Element.h"
#include "flex_mesh/mesh/Mesh.h"
#include "flex_mesh/maps/AdjacencyLists.h"

#include "second_party_metis/metis.h"

using namespace flex_mesh;

int main()
{
	initialize_flex_mesh(1e-14);

	using NodeTyped = Node;
	using ElementTyped = Element<NodeTyped>;
	using MeshType = Mesh<NodeTyped, ElementTyped>;

	std::vector<double> coordinates(10);
	std::iota(coordinates.begin(), coordinates.end(), 0.0);

	MeshType mesh = flex_mesh::make_grid<MeshType>(
		coordinates, coordinates, coordinates, false);

	flex_mesh::partition_metis(mesh, 3);
	flex_mesh::write_vtk(mesh, "total_mesh.vtk", true, false);

	auto mesh_part0 = flex_mesh::extract_partition_mesh(mesh, 0);
	flex_mesh::write_vtk(mesh_part0, "partition_0.vtk", true, false);

	return 0;
}