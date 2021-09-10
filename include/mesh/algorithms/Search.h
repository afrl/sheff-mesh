#pragma once

#include "mesh/Mesh.h"
#include "mesh/surface/Surface.h"

#include <utility>

namespace mesh
{
/// @brief Finds extent of a list of points (most negative and positive point forming the bounding
/// box of the mesh).
/// @returns A pair of points with most negative point being first and most positive being second.
std::pair<Point, Point> find_extent(const std::vector<simple_math::Point>& points);

/// @brief Finds extent of the mesh (most negative and positive point forming the bounding box of
/// the mesh).
/// @returns A pair of points with most negative point being first and most positive being second.
std::pair<Point, Point> find_extent(const Mesh& mesh);

/// @brief Returns the extent of multiple extents.
std::pair<Point, Point> find_extent(const std::vector<std::pair<Point, Point>>& extents);

/// @brief Finds all elements whose centroid is the in the bounding box defined by (min, max)
std::vector<LocalIndex> find_elements_in_box(const Mesh& mesh, const Point& min, const Point& max);

/// @brief Finds global extent of the mesh (most negative and positive point forming the bounding
/// box of the mesh) taking into considering that the mesh may be distributed across ranks.
/// @returns A pair of points with most negative point being first and most positive being second.
std::pair<Point, Point> find_global_extent(const Mesh& mesh);

/// @brief Detects if two bounding boxes that are aligned with the coordinate axes intersect
///
/// The algorithm leverages the separate axis theorem.
bool do_bounding_boxes_intersect(const std::pair<Point, Point>& bb1,
                                 const std::pair<Point, Point>& bb2);

/// @brief Detects if two bounding boxes that are aligned with the coordinate axes intersect
///
/// The algorithm leverages the separate axis theorem.
bool do_bounding_boxes_intersect(const std::array<Point, 2>& bb1,
                                 const std::array<Point, 2>& bb2);

double bounding_box_volume(const std::pair<Point, Point>& bb);
double bounding_box_volume(const std::array<Point, 2>& bb);

/// @brief Finds the maximum component ID in the mesh (ignores other ranks)
ComponentIndex find_max_component_id(const Mesh& mesh);

/// @brief Finds the maximum component ID in the mesh taking into considering that the mesh may be
/// distributed across ranks.
ComponentIndex find_global_max_component_id(const Mesh& mesh);

struct BoundaryInfo
{
    MeshSurface surface;
    std::vector<LocalIndex> parent_volume_elements;
    std::vector<ConnIndex> face_in_parent_volume_elements;
};

/// @brief TODO write doc string
/// @param mesh
BoundaryInfo extract_boundary(const Mesh& mesh);

/// @brief TODO write doc string
/// @param mesh
/// @param elements_to_search
BoundaryInfo extract_boundary(const Mesh& mesh, const std::vector<LocalIndex>& elements_to_search);

/// @brief Extracts a submesh for a single component from a mesh that may contain multiple
/// components.
///
/// Note that this is different from extract_boundary, which extracts the boundary a given set of
/// elements.  This simply returns a mesh containing the elements that have the corresponding
/// component ID (and only the nodes connected to the respective elements).
Mesh extract_component(const Mesh& mesh, const ComponentIndex& component_id);

/// @brief Simple struct to hold result of a call to find_unmatched_edges
/// @see find_unmatched_edges
struct UnmatchedEdgeInfo
{
    std::vector<LocalIndex> parent_element_indices;
    std::vector<EdgeType> edge_types;
    std::vector<std::vector<LocalIndex>> edge_connectivities;
};

/// Finds all edges that appear in a single element and returns information about them.
UnmatchedEdgeInfo find_unmatched_edges(const Mesh& mesh);

/// A helper struct for find_faces_from_node_set
struct ElementFacePair
{
    LocalIndex local_element_index;
    ConnIndex face_index;
};

/// @brief Find all faces such that all nodes in the face appears in the given node set.
std::vector<ElementFacePair> find_faces_from_node_set(const Mesh& mesh,
                                                      const std::set<LocalIndex>& node_set);

/// @brief Find an element that contains the given point
/// @returns A pair indicating (first) if an element was found and (second) the element index
std::pair<bool, LocalIndex> find_element_with_point(const Mesh& mesh, const simple_math::Point x);

/// @brief Counts how many elements of each unique element type in a mesh.
std::map<ElementType, LocalIndex> count_element_types(const Mesh& mesh);
}