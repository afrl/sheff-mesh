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
	using FaceTyped = Face<DistribNode<Node>>;

	NodeTyped na(0, 0, { 0, 0, 0 });
	NodeTyped nb(1, 1, { 1, 0, 0 });
	NodeTyped nc(2, 2, { 1, 1, 0 });
	NodeTyped nd(2, 2, { 0, 1, 0 });

	FaceTyped a(FaceType::Q3D4N, { &na, &nb, &nc, &nd });
	assert(a.get_type() == FaceType::Q3D4N);
	for (auto n_ptr : a)
		assert(n_ptr->get_parent_partition() == 0);
	// Technically, this could be an invalid assert, since the hash of the node addresses
	// actually be exactly 0, but the changes are astronomically small.
	assert(a.get_conn_hash_id() > 0);

	Point mean = a.get_mean_coordinate();
	assert(distance(mean, { 0.5, 0.5, 0 }) < 1e-16);

	auto edges = a.get_edges();
	assert(edges.size() == 4);

	auto outward_norm = a.get_outward_normal();
	assert((outward_norm - Vector3D(0, 0, 1)).norm() < 1e-16);
	double projection = a.calc_signed_norm_projection({ 0.5, 0.5, 1 });
	assert(fabs(projection - 1.0) < 1e-16);

	return 0;
}