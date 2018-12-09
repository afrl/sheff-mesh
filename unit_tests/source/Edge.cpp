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
	using EdgeTyped = Edge<DistribNode<Node>>;

	NodeTyped na(0, 0, { 0, 0, 0 });
	NodeTyped nb(1, 1, { 0.5, 0, 0 });
	EdgeTyped a(EdgeType::E2N, { &na, &nb });
	assert(a.get_type() == EdgeType::E2N);
	for (auto n_ptr : a)
		assert(n_ptr->get_parent_partition() == 0);
	// Technically, this could be an invalid assert, since the hash of the node addresses
	// actually be exactly 0, but the changes are astronomically small.
	assert(a.get_conn_hash_id() > 0);
	assert(fabs(a.get_length() - 0.5) < 1e-16);
	Point mean = a.get_mean_coordinate();
	assert(distance(mean, { 0.25, 0, 0 }) < 1e-16);

	NodeTyped nc(2, 2, { 1, 0, 0 });
	a.swap_nodes(&nb, &nc);
	assert(fabs(a.get_length() - 1) < 1e-16);

	EdgeTyped b(EdgeType::E2N, { &na, &nb });
	assert(b < a);

	return 0;
}