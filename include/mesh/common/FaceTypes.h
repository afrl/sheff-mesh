#pragma once

#include "EdgeTypes.h"
#include "mesh/Defines.h"
#include "utility/Exceptions.h"

#include <vector>

namespace mesh
{
/// Describes the different types of faces in terms of the geometry/topology.
/// Note: faces are a bit weird in 1D.  In 1D, a face will consist of a single node.
enum class FaceType : char
{
    UNKNGEOM = 0, ///< Used if a conversion fails
    // 1D types 10-19
    L1D1N = 11, ///< 1-D Linear 1-node
    // 2D types 20-24 is linear and 25-29 is quadratic
    L2D2N = 21, ///< 2-D 2-node face (linear)
    L2D3N = 22, ///< 2-D 3-node face (quadratic)
    // Linear 3D 30-39
    T3D3N = 31, ///< triangular 3-node face (linear)
    Q3D4N = 32, ///< quadrilateral 4-node face (linear)
    // Quadratic 3D 41-49
    T3D6N = 41, ///< triangular 6-node face (quadratic)
    Q3D8N = 42, ///< quadrilateral 8-node face (quadratic)
};

/// Returns the number of nodes that the geometry type will have.
static ConnIndex get_num_nodes_for_type(const FaceType& geom_type);
/// Given an face geometry type, this function will fill the vectors:
/// - types: with the types of edges
/// - conn: with connectivity of each face relative to the face connectivity
/// \attention Will overwrite anything in the vectors!
static const std::vector<EdgeInfo> get_edge_info(const FaceType& geom_type);
/// Given a geometry type, this function will return the linear version
/// if one exists.  Otherwise it will return UNKNGEOM.
static FaceType get_linearized_type(const FaceType& geom_type);
/// Returns the permuation vector to go from the quadratic connectivity
/// to the linear connectivity.  If the type is already linear, then
/// an empty vector is returned.  Will throw for unrecognized types.
//static std::vector<ConnIndex> get_linearized_permutation(const FaceType& geom_type);
/// Returns a bool indicating if the geometry type is linear or not
static bool is_type_linear(const FaceType& geom_type);
/// Returns the permutation for flipping an element's connectivity (reversing outward normal).
static std::vector<ConnIndex> flip_conn_permute(const FaceType& geom_type);


/// @brief A simple struct to hold information about a face that is part of a element.
struct FaceInfo
{
    FaceType type;
    std::vector<ConnIndex> edges;
    std::vector<ConnIndex> connectivity;
};

/// @brief A simple struct to hold information about a describes a face as a polygonal facet.
/// Returned by mesh::get_facet_info.
struct FacetInfo
{
    /// Clock-wise points of the facet (relative to element connectivity)
    std::vector<ConnIndex> points;
    /// The edge index that the segment is a part of
    std::vector<int8_t> segment_markers;
};

//////////////////
// Definitiions //
//////////////////

ConnIndex get_num_nodes_for_type(const FaceType& geom_type)
{
    if (geom_type == FaceType::L1D1N)
        return 1;
    else if (geom_type == FaceType::L2D2N)
        return 2;
    else if (geom_type == FaceType::L2D3N)
        return 3;
    else if (geom_type == FaceType::T3D3N)
        return 3;
    else if (geom_type == FaceType::Q3D4N)
        return 4;
    else if (geom_type == FaceType::T3D6N)
        return 6;
    else if (geom_type == FaceType::Q3D8N)
        return 8;
    else
        throw _UTIL_MSG_EXCEPT("Face type not supported");
}

const std::vector<EdgeInfo> get_edge_info(const FaceType& geom_type)
{
    switch (geom_type)
    {
        // 2D Connectivities
        case FaceType::L2D2N: return { { EdgeType::E2N, { 0, 1 } } };
        case FaceType::L2D3N: return { { EdgeType::E3N, { 0, 1, 2 } } };
        // 3D Connectivities
        case FaceType::T3D3N:
            return { { EdgeType::E2N, { 0, 1 } },
                     { EdgeType::E2N, { 1, 2 } },
                     { EdgeType::E2N, { 2, 0 } } };
        case FaceType::Q3D4N:
            return { { EdgeType::E2N, { 0, 1 } },
                     { EdgeType::E2N, { 1, 2 } },
                     { EdgeType::E2N, { 2, 3 } },
                     { EdgeType::E2N, { 3, 0 } } };
        case FaceType::T3D6N:
            return { { EdgeType::E3N, { 0, 1, 2 } },
                     { EdgeType::E3N, { 2, 3, 4 } },
                     { EdgeType::E3N, { 4, 5, 0 } } };
        case FaceType::Q3D8N:
            return { { EdgeType::E3N, { 0, 1, 2 } },
                     { EdgeType::E3N, { 2, 3, 4 } },
                     { EdgeType::E3N, { 4, 5, 6 } },
                     { EdgeType::E3N, { 6, 7, 0 } } };
        default: throw _UTIL_MSG_EXCEPT("Edge connectivity not implemented.");
    }
}

FaceType get_linearized_type(const FaceType& geom_type)
{
    // 2D
    if (geom_type == FaceType::L2D3N)
        return FaceType::L2D2N;
    // 3D
    else if (geom_type == FaceType::Q3D8N)
        return FaceType::Q3D4N;
    else if (geom_type == FaceType::T3D6N)
        return FaceType::T3D3N;
    else
        return FaceType::UNKNGEOM;
}

bool is_type_linear(const FaceType& geom_type)
{
    return geom_type == FaceType::L1D1N || geom_type == FaceType::L2D2N ||
           geom_type == FaceType::Q3D4N || geom_type == FaceType::T3D3N;
}

std::vector<ConnIndex> flip_conn_permute(const FaceType& geom_type)
{
    switch (geom_type)
    {
        // 2D faces
        case FaceType::L2D2N: return { 1, 0 };
        case FaceType::L2D3N: return { 2, 1, 0 };
        // 3D faces
        case FaceType::T3D3N: return { 2, 1, 0 };
        case FaceType::Q3D4N: return { 3, 2, 1, 0 };
        case FaceType::T3D6N: return { 0, 5, 4, 3, 2, 1 };
        case FaceType::Q3D8N: return { 0, 7, 6, 5, 4, 3, 2, 1 };
        default: throw _UTIL_MSG_EXCEPT("Flip is not implemented for face type.");
    }
}
};
