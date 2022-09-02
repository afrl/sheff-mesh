#include "mesh/algorithms/Search.h"

#include "mesh/basis/BasisFunctions.h"
#include "mesh/face/Face.h"
#include "mesh/surface/Surface.h"
#include "simple_math/Utility.h"
#include "utility/Hash.h"
#include "utility/MPIUtils.h"

#include <unordered_map>


namespace mesh
{
std::pair<Point, Point> find_extent(const std::vector<simple_math::Point>& points)
{
    Point neg_point = { INFINITY, INFINITY, INFINITY };
    Point pos_point = { -INFINITY, -INFINITY, -INFINITY };
    for (auto& p : points)
    {
        for (auto i = 0; i < neg_point.size(); i++)
        {
            if (p[i] < neg_point[i])
            {
                neg_point[i] = p[i];
            }
            if (p[i] > pos_point[i])
            {
                pos_point[i] = p[i];
            }
        }
    }
    return { neg_point, pos_point };
}

std::pair<Point, Point> find_extent(const Mesh& mesh)
{
    Point neg_point = { INFINITY, INFINITY, INFINITY };
    Point pos_point = { -INFINITY, -INFINITY, -INFINITY };
    for (auto node : mesh.nodes())
    {
        auto p = node.point();
        for (auto i = 0; i < neg_point.size(); i++)
        {
            if (p[i] < neg_point[i])
            {
                neg_point[i] = p[i];
            }
            if (p[i] > pos_point[i])
            {
                pos_point[i] = p[i];
            }
        }
    }
    return { neg_point, pos_point };
}
std::pair<Point, Point> find_extent(const std::vector<std::pair<Point, Point>>& extents)
{
    Point neg_point = { INFINITY, INFINITY, INFINITY };
    Point pos_point = { -INFINITY, -INFINITY, -INFINITY };
    for (auto& extent : extents)
    {
        for (auto& p : { extent.first, extent.second })
        {
            for (auto i = 0; i < neg_point.size(); i++)
            {
                if (p[i] < neg_point[i])
                {
                    neg_point[i] = p[i];
                }
                if (p[i] > pos_point[i])
                {
                    pos_point[i] = p[i];
                }
            }
        }
    }
    return { neg_point, pos_point };
}

std::vector<LocalIndex> find_nodes_in_box(const Mesh& mesh, const Point& min, const Point& max)
{
    std::vector<LocalIndex> node_indices;
    for (auto n : mesh.nodes())
    {
        bool inside = true;
        for (auto i : utility::range(3))
        {
            inside &= n.point(i) >= min[i] && n.point(i) <= max[i];
        }
        if (inside)
        {
            node_indices.push_back(n.local_index());
        }
    }
    return node_indices;
}

std::vector<LocalIndex> find_elements_in_box(const Mesh& mesh, const Point& min, const Point& max)
{
    std::vector<LocalIndex> element_indices;
    for (auto e : mesh.elements())
    {
        auto m      = calc_mean_coordinate(e);
        bool inside = true;
        for (auto i : utility::range(3))
        {
            inside &= m[i] >= min[i] && m[i] <= max[i];
        }
        if (inside)
        {
            element_indices.push_back(e.local_index());
        }
    }
    return element_indices;
}

std::pair<Point, Point> find_global_extent(const Mesh& mesh)
{
    auto local_extent = find_extent(mesh);
    std::pair<Point, Point> global_extent;
    for (auto i : utility::range(3))
    {
        global_extent.first[i]  = utility::MPIUtils::min_across_ranks(local_extent.first[i]);
        global_extent.second[i] = utility::MPIUtils::max_across_ranks(local_extent.second[i]);
    }
    return global_extent;
}

namespace detail
{
/// @brief Returns whether line segments a and b overlap in anyway
/// Detects the following cases:
///
/// a: x-----x
/// b:    x-----x
///
/// a:     x-----x
/// b: x-----x
///
/// a: x-----------x
/// b:    x-----x
///
/// a:    x-----x
/// b: x-----------x
bool do_segments_overlap(const double& a1, const double& a2, const double& b1, const double& b2)
{
    const double t = 1e-6 * std::min(a2 - a1, b2 - b1);
    if ((a1 >= b1 - t && a1 <= b2 + t) || (a2 >= b1 - t && a2 <= b2 + t) ||
        (b1 >= a1 - t && b1 <= a2 + t) || (b2 >= a1 - t && b2 <= a2 + t))
    {
        return true;
    }
    else
    {
        return false;
    }
}
}

bool do_bounding_boxes_intersect(const std::pair<Point, Point>& bb1,
                                 const std::pair<Point, Point>& bb2)
{
    if (detail::do_segments_overlap(bb1.first[0], bb1.second[0], bb2.first[0], bb2.second[0]) &&
        detail::do_segments_overlap(bb1.first[1], bb1.second[1], bb2.first[1], bb2.second[1]) &&
        detail::do_segments_overlap(bb1.first[2], bb1.second[2], bb2.first[2], bb2.second[2]))
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool do_bounding_boxes_intersect(const std::array<Point, 2>& bb1, const std::array<Point, 2>& bb2)
{
    if (detail::do_segments_overlap(bb1[0][0], bb1[1][0], bb2[0][0], bb2[1][0]) &&
        detail::do_segments_overlap(bb1[0][1], bb1[1][1], bb2[0][1], bb2[1][1]) &&
        detail::do_segments_overlap(bb1[0][2], bb1[1][2], bb2[0][2], bb2[1][2]))
    {
        return true;
    }
    else
    {
        return false;
    }
}

double bounding_box_volume(const std::pair<Point, Point>& bb)
{
    auto bb_p1 = simple_math::point_as_vector(bb.first);
    auto bb_p2 = simple_math::point_as_vector(bb.second);
    return (bb_p2 - bb_p1).norm();
}

double bounding_box_volume(const std::array<Point, 2>& bb)
{
    auto bb_p1 = simple_math::point_as_vector(bb[0]);
    auto bb_p2 = simple_math::point_as_vector(bb[1]);
    return (bb_p2 - bb_p1).norm();
}

ComponentIndex find_max_component_id(const Mesh& mesh)
{
    mesh::ComponentIndex curr_max_id = 0;
    for (auto& element : mesh.elements())
    {
        if (element.component_id() > curr_max_id)
        {
            curr_max_id = element.component_id();
        }
    }
    return curr_max_id;
}

ComponentIndex find_global_max_component_id(const Mesh& mesh)
{
    return utility::MPIUtils::max_across_ranks(find_max_component_id(mesh));
}

BoundaryInfo extract_boundary(const Mesh& mesh)
{
    std::vector<LocalIndex> elements_to_search(mesh.elements().size());
    std::iota(elements_to_search.begin(), elements_to_search.end(), 0);
    return extract_boundary(mesh, elements_to_search);
}

BoundaryInfo extract_boundary(const Mesh& mesh, const std::vector<LocalIndex>& elements_to_search)
{
    // A map from the face hash to the pair of (element index, face index) for all unmatched faces.
    // If a match is found, then the existing object is removed (thus assuming no more than 2
    // elements share a face).
    std::unordered_map<std::size_t, std::pair<LocalIndex, ConnIndex>> unmatched_faces;
    for (auto& element_index : elements_to_search)
    {
        FaceContainer faces(mesh.elements()[element_index]);
        for (auto& face : faces)
        {
            auto hash      = utility::hash_sorted(face.connectivity());
            auto find_iter = unmatched_faces.find(hash);
            if (find_iter == unmatched_faces.end())
            {
                unmatched_faces[hash] = { element_index, face.index() };
            }
            else
            {
                unmatched_faces.erase(find_iter);
            }
        }
    }

    BoundaryInfo boundary_info = { MeshSurface(mesh), {} };
    for (auto& p : unmatched_faces)
    {
        auto& element_index = p.second.first;
        auto& face_index    = p.second.second;

        FaceContainer faces(mesh.elements()[element_index]);
        auto face = faces[face_index];
        boundary_info.surface.add(face);
        boundary_info.parent_volume_elements.push_back(element_index);
        boundary_info.face_in_parent_volume_elements.push_back(face_index);
    }
    return boundary_info;
}

Mesh extract_component(const Mesh& mesh, const ComponentIndex& component_id)
{
    std::vector<LocalIndex> elements_in_component;
    for (auto e : mesh.elements())
    {
        elements_in_component.push_back(e.local_index());
    }

    Mesh component_mesh;
    component_mesh.copy(mesh, elements_in_component);
    return component_mesh;
}

UnmatchedEdgeInfo find_unmatched_edges(const Mesh& mesh)
{
    // Will store the connectivity of an edge in terms of the node numbering in the mesh
    std::vector<LocalIndex> edge_conn;
    // A map from the hash of an edge's sorted connectivity to the number of elements that contain
    // the edge, the first element that contained it (based on local element numbering), the type of
    // edge, and the connectivity of the edge.
    std::unordered_map<std::size_t, std::tuple<int, LocalIndex, EdgeType, std::vector<LocalIndex>>>
      edge_hash_to_count_and_info;
    for (auto e : mesh.elements())
    {
        auto e_conn = e.connectivity();
        for (auto& edge : get_edge_info(e.type()))
        {
            edge_conn.resize(edge.connectivity.size());
            for (auto i : utility::indices(edge.connectivity))
            {
                edge_conn[i] = e_conn[edge.connectivity[i]];
            }
            auto edge_hash = utility::hash_sorted(edge_conn);
            auto find_iter = edge_hash_to_count_and_info.find(edge_hash);
            if (find_iter == edge_hash_to_count_and_info.end())
            {
                edge_hash_to_count_and_info[edge_hash] = { 1,
                                                           e.local_index(),
                                                           edge.type,
                                                           edge_conn };
            }
            else
            {
                std::get<0>(find_iter->second)++;
            }
        }
    }

    UnmatchedEdgeInfo unmatched_info;
    for (auto& key_pair : edge_hash_to_count_and_info)
    {
        // The unmatched edges only appear in one element
        if (std::get<0>(key_pair.second) == 1)
        {
            unmatched_info.parent_element_indices.push_back(std::get<1>(key_pair.second));
            unmatched_info.edge_types.push_back(std::get<2>(key_pair.second));
            unmatched_info.edge_connectivities.push_back(std::get<3>(key_pair.second));
        }
    }
    return unmatched_info;
}

std::vector<ElementFacePair> find_faces_from_node_set(const Mesh& mesh,
                                                      const std::set<LocalIndex>& node_set)
{
    std::vector<ElementFacePair> faces_in_node_set;
    for (auto& element : mesh.elements())
    {
        FaceContainer faces(element);
        for (auto face_index : utility::indices(faces))
        {
            auto face_conn = faces[face_index].connectivity();
            if (std::all_of(face_conn.begin(),
                            face_conn.end(),
                            [&node_set](const LocalIndex& i) -> bool
                            { return node_set.count(i) != 0; }))
            {
                faces_in_node_set.push_back({ element.local_index(), face_index });
            }
        }
    }
    return faces_in_node_set;
}

std::pair<bool, LocalIndex> find_element_with_point(const Mesh& mesh, const simple_math::Point x)
{
    for (auto& element : mesh.elements())
    {
        if (contains_point(element, x))
        {
            return { true, element.local_index() };
        }
    }
    return { false, 0 };
}

std::map<ElementType, LocalIndex> count_element_types(const Mesh& mesh)
{
    std::map<ElementType, LocalIndex> type_count_map;
    for (auto& element : mesh.elements())
    {
        auto find_iter = type_count_map.find(element.type());
        if (find_iter == type_count_map.end())
        {
            type_count_map[element.type()] = 1;
        }
        else
        {
            find_iter->second++;
        }
    }
    return type_count_map;
}

simple_math::Point find_parametric_point(const ElementConstAccessor& element,
                                         const simple_math::Point& global_point)
{
    auto global_coord = simple_math::point_as_vector(global_point);
    simple_math::GeneralMatrix node_coords(element.connectivity().size(), 3);
    for (auto i : utility::range(element.connectivity().size()))
    {
        node_coords.row(i) = element.node(i).coord();
    }

    simple_math::GeneralVector values;
    simple_math::GeneralMatrix derivatives;
    simple_math::Vector3D curr_guess = simple_math::Vector3D::Zero();
    simple_math::Point curr_guess_p;
    int curr_iteration = 0;

    do
    {
        curr_guess_p = simple_math::to_point(curr_guess);
        BasisFunctions::values(element.type(), curr_guess_p, values);
        auto residual = node_coords.transpose() * values - global_coord;
        if (residual.norm() < 1e-8)
        {
            return curr_guess_p;
        }

        BasisFunctions::derivatives(element.type(), curr_guess_p, derivatives);
        auto jacobian = derivatives.transpose() * node_coords;

        curr_guess -= jacobian.inverse() * residual;
        curr_iteration++;
    } while (curr_iteration < 20);

    throw _UTIL_MSG_EXCEPT("Inverse problem to find parametric coordinate did not converge.");
}
}