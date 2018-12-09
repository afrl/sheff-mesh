#include <chrono>
#include <functional>
#include <iostream>
#include <iomanip>

#include "flex_mesh/common/Common.h"
#include "flex_mesh/common/Connectivity.h"
#include "flex_mesh/node/Node.h"
#include "flex_mesh/node/DistribNode.h"
#include "flex_mesh/node/NodeUtilities.h"
#include "flex_mesh/edge/Edge.h"
#include "flex_mesh/face/Face.h"
#include "flex_mesh/element/Element.h"
#include "flex_mesh/element/DistribElement.h"
#include "flex_mesh/element/MeshGenElement.h"

using namespace flex_mesh;

int main()
{
	initialize_flex_mesh(1e-14);

	using NodeTyped = DistribNode<Node>;
	using ElementTyped = MeshGenElement<NodeTyped, DistribElement<NodeTyped, Element<NodeTyped>>>;
	using MeshType = Mesh<NodeTyped, ElementTyped>;

	MeshType mesh;
	NodeTyped * na = mesh.add_node_by_coord({ 0, 0, 0 });
	NodeTyped * nb = mesh.add_node_by_coord({ 1, 0, 0 });
	NodeTyped * nc = mesh.add_node_by_coord({ 1, 1, 0 });
	NodeTyped * nd = mesh.add_node_by_coord({ 0, 1, 0 });
	ElementTyped * ea = mesh.add_element();
	ea->initialize(ElementType::Q3D4N, { na, nb, nc, nd });
	mesh.update_global_numbering();
	mesh.sort_by_partition();

	// TODO: add asserts

	return 0;
}