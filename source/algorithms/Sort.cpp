#include "mesh/algorithms/Sort.h"

#include "simple_math/Utility.h"
#include "utility/Permutation.h"

namespace mesh
{
void sort_node_indices_by_coord(const Mesh& mesh, std::vector<mesh::LocalIndex>& node_indices)
{
    auto cmp = [&mesh](const mesh::LocalIndex& a, const mesh::LocalIndex& b)
    {
        const auto& A = mesh.nodes()[a];
        const auto& B = mesh.nodes()[b];
        return (A.x() < B.x()) || ((A.y() < B.y()) && (A.x() == B.x())) ||
               ((A.z() < B.z()) && (A.x() == B.x()) && (A.y() == B.y()));
    };
    std::sort(node_indices.begin(), node_indices.end(), cmp);
}

void sort_element_face_pairs_by_coord(const Mesh& mesh,
                                      std::vector<ElementFacePair>& element_face_pairs)
{
    using namespace simple_math::point_operators;

    std::vector<Point> face_mean_coords;
    face_mean_coords.reserve(element_face_pairs.size());
    for (auto& ef_pair : element_face_pairs)
    {
        const auto& element = mesh.elements()[ef_pair.local_element_index];
        FaceContainer faces(element);
        face_mean_coords.push_back(calc_mean_coordinate(faces[ef_pair.face_index]));
    }

    auto cmp = [&mesh](const Point& A, const Point& B)
    {
        return (A[0] < B[0]) || ((A[1] < B[1]) && (A[0] == B[0])) ||
               ((A[2] < B[2]) && (A[0] == B[0]) && (A[1] == B[1]));
    };
    auto backward_map = utility::find_sorted_permutation(face_mean_coords, cmp);
    utility::apply_permutation(element_face_pairs, backward_map);
}
}