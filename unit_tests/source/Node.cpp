#include <chrono>
#include <functional>
#include <iostream>
#include <iomanip>

#include "flex_mesh/common/Common.h"
#include "flex_mesh/common/Connectivity.h"
#include "flex_mesh/node/Node.h"
#include "flex_mesh/node/NodeUtilities.h"

using namespace flex_mesh;

int main()
{
	initialize_flex_mesh(1e-14);

	using NodeTyped = Node;

	NodeTyped a({ 1, 0, 0 });
	assert(a.get_global_number() == 0);
	assert(a.get_coord() == Point({ 1, 0, 0 }));
	assert(a.get_parent_partition() == 0);

	NodeTyped b(1, { 0, 0, 1 }, 1);
	assert(b.get_global_number() == 1);
	assert(b.get_coord() == Point({ 0, 0, 1 }));
	assert(b.get_parent_partition() == 1);

	assert(a > b);
	assert(b < a);
	assert(!(a == b));
	assert(NodeTyped::num_based_compare(a, b));
	assert(NodeTyped::part_based_compare(a, b));

	Node* test1 = &a;
	Node** test2 = new Node*;
	*test2 = &a;

	std::array<NodeTyped *, 2> nodes = { &a, &b };
	Point mean = get_mean_coordinate_of_nodes(nodes);
	assert(distance(mean, { 0.5, 0, 0.5 }) < 1e-16);
	assert(fabs(distance(a, b) - sqrt(2)) < 1e-16);

	return 0;
}