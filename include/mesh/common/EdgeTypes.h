#pragma once

#include "utility/Exceptions.h"
#include "mesh/Defines.h"

#include <vector>

namespace mesh
{
	/// Describes the different types of edges in terms of the geometry/topology.
	/// Note: edges are a bit weird in 2D or 1D.  In 1D/2D the edge will be same
	///       as a 1D/2D face basically.
	enum class EdgeType : char {
		UNKNGEOM = 0, ///< Used if a conversion fails
		E2N = 10, ///< 2-node edge
		E3N = 11, ///< 3-node edge
	};

	/// Returns the number of nodes that the geometry type will have.
	static ConnIndex get_num_nodes_for_type(const EdgeType& geom_type);

    /// A simple struct to hold information about an edge that is part of a face/element.
    struct EdgeInfo
    {
        EdgeType type;
        std::vector<ConnIndex> connectivity;
    };

	//////////////////
	// Definitiions //
	//////////////////

    ConnIndex get_num_nodes_for_type(const EdgeType & geom_type)
	{
		if (geom_type == EdgeType::E2N) return 2;
		else if (geom_type == EdgeType::E3N) return 3;
		else throw _UTIL_MSG_EXCEPT("Edge type not supported");
	}
};
