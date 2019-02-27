#include <chrono>
#include <functional>
#include <iostream>
#include <iomanip>
#include <memory>

#define EIGEN_MAX_CPP_VER 11

#include "flex_mesh/common/Common.h"
#include "flex_mesh/common/Connectivity.h"
#include "flex_mesh/node/Node.h"
#include "flex_mesh/node/NodeUtilities.h"

using namespace flex_mesh;

int main()
{
	initialize_flex_mesh(1e-14);

	using NodeTyped = Node;

	NodeTyped na(0, { 0, 0, 0 });
	NodeTyped nb(1, { 0.5, 0, 0 });
	std::vector<NodeTyped*> c_vec;

	{
		Connectivity<NodeTyped, uint16_t> c({ &na, &nb });
		assert(*c.begin() == &na);
		assert(c.size() == 2);
		for (auto n_ptr : c)
			assert(n_ptr->get_parent_partition() == 0);
		assert(c[0]->get_global_number() == 0);
		assert(c[1]->get_global_number() == 1);

		c_vec = c.copy_as_vector();
		assert(c_vec[0] == c[0]);
		assert(&c_vec[0] != c.begin());
	}

	assert(c_vec[0]->get_global_number() == 0);
	assert(c_vec[1]->get_global_number() == 1);

	return 0;
}
