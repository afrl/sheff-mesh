#include "mesh/algorithms/Search.h"

#include "mesh/face/Face.h"
#include "mesh/surface/Surface.h"
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
    for (auto& node : mesh.nodes())
    {
        auto& p = node.point();
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

ComponentIndex find_global_max_component_id(const Mesh& mesh)
{
    mesh::ComponentIndex curr_max_id = 0;
    for (auto& element : mesh.elements())
    {
        if (element.component_id() > curr_max_id)
        {
            curr_max_id = element.component_id();
        }
    }
    return utility::MPIUtils::max_across_ranks(curr_max_id);
}

MeshSurface extract_boundary(Mesh& mesh)
{
    std::vector<LocalIndex> elements_to_search(mesh.elements().size());
    std::iota(elements_to_search.begin(), elements_to_search.end(), 0);
    return extract_boundary(mesh, elements_to_search);
}

MeshSurface extract_boundary(Mesh& mesh, const std::vector<LocalIndex>& elements_to_search)
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

    MeshSurface surface_elements(mesh);
    for (auto& p : unmatched_faces)
    {
        auto& element_index = p.second.first;
        auto& face_index    = p.second.second;

        FaceContainer faces(mesh.elements()[element_index]);
        auto& face = faces[face_index];
        surface_elements.add(face);
    }
    return surface_elements;
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
            auto& face_conn = faces[face_index].connectivity();
            if (std::all_of(
                  face_conn.begin(),
                  face_conn.end(),
                  [&node_set](const LocalIndex& i) -> bool { return node_set.count(i) != 0; }))
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
}