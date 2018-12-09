#pragma once

#include "flex_mesh/node/Node.h"

namespace flex_mesh
{
	/// Sets all the static variables in the library.  This should be called 
	/// once (and only once) by the program using this library.
	static void initialize_flex_mesh(
		const double node_compare_tol)
	{
		Node::set_compare_tolerance(node_compare_tol);
	}
};