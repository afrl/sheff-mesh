#include "mesh/algorithms/Search.h"

#include "mesh/face/Face.h"
#include "mesh/surface/Surface.h"
#include "utility/Hash.h"

#include <unordered_map>

namespace mesh
{
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
}