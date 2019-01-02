#include <chrono>
#include <functional>
#include <iostream>
#include <iomanip>

#include "flex_mesh/common/Common.h"
#include "flex_mesh/common/Connectivity.h"
#include "flex_mesh/node/Node.h"
#include "flex_mesh/node/NodeUtilities.h"
#include "flex_mesh/element/Element.h"
#include "flex_mesh/element/MeshGenElement.h"

using namespace flex_mesh;

int main()
{
	initialize_flex_mesh(1e-14);

	using NodeTyped = Node;
	using ElementTyped = MeshGenElement<Element<NodeTyped>>;

	NodeTyped na(0, { 0, 0, 0 });
	NodeTyped nb(1, { 1, 0, 0 });
	NodeTyped nc(2, { 1, 1, 0 });
	NodeTyped nd(2, { 0, 1, 0 });

	ElementTyped a;
	a.initialize(ElementType::Q2D4N, { &na, &nb, &nc, &nd });
	a.set_component_id(1);
	assert(a.get_component_id() == 1);
	a.set_global_number(1);
	assert(a.get_global_number() == 1);
	a.set_parent_partition(1);
	assert(a.get_parent_partition() == 1);

	assert(a.contains_node(na));
	assert(a.get_conn_hash_id() > 0);
	Point mean = a.get_mean_coordinate();
	assert(distance(mean, { 0.5, 0.5, 0 }) < 1e-16);

	assert(a.get_lcs().size() == 0);
	a.allocate_lcs(LCSAssignment::None);
	assert(a.get_lcs().size() == 0);
	a.allocate_lcs(LCSAssignment::Element);
	assert(a.get_lcs().size() == 1);
	a.allocate_lcs(LCSAssignment::Component);
	assert(a.get_lcs().size() == 1);
	a.allocate_lcs(LCSAssignment::Node);
	assert(a.get_lcs().size() == 4);

	ElementVector S;
	ElementMatrix3D dS;
	a.get_shape_function_values({ 0, 0, 0 }, S);
	a.get_shape_function_derivatives({ 0, 0, 0 }, dS);

	return 0;
}