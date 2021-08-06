#include "mesh/algorithms/Refinement.h"

#include "mesh/algorithms/Manipulation.h"
#include "simple_math/Utility.h"
#include "utility/Hash.h"

using namespace simple_math::point_operators;

namespace mesh
{
void refine_linear_tetrahedra(Mesh& mesh,
                              const ElementAccessor& element,
                              std::map<std::size_t, LocalIndex>& edge_hash_to_new_node_map)
{
    // Based on connectivity: (with new nodes labeled)
    //      3
    //     /|\
    //    / | \
    //  7/  |8 \9
    //  /   |   \
    // 0----1----2
    //   4  6  5
    //
    // View of bottom face:
    //      6
    // 0---------2
    //  \       /
    //   \     /
    //   4\   /5
    //     \ /
    //      1
    

    auto orig_conn = element.connectivity();
    auto edge_info = get_edge_info(element.type());

    // Find or create new midside nodes (4-9 in the figure above)
    std::vector<mesh::LocalIndex> c =
      std::vector<mesh::LocalIndex>(orig_conn.begin(), orig_conn.end());
    for (auto& edge : edge_info)
    {
        assert(edge.type == EdgeType::E2N);

        std::vector<LocalIndex> edge_conn;
        edge_conn.reserve(edge.connectivity.size());
        for (auto& c : edge.connectivity)
        {
            edge_conn.push_back(element.connectivity()[c]);
        }
        auto edge_hash = utility::hash_sorted(edge_conn);

        auto find_iter = edge_hash_to_new_node_map.find(edge_hash);
        if (find_iter == edge_hash_to_new_node_map.end())
        { // If it doesn't exist, add a midside node
            auto new_coord =
              (mesh.nodes()[edge_conn[1]].coord() + mesh.nodes()[edge_conn[0]].coord()) / 2.0;
            mesh.nodes().add(simple_math::to_point(new_coord));
            edge_hash_to_new_node_map[edge_hash] = mesh.nodes().size() - 1;
            c.push_back(mesh.nodes().size() - 1);
        }
        else
        {
            c.push_back(find_iter->second);
        }
    }

    // Create the new tetrahedra
    // See https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.308.2549&rep=rep1&type=pdf
    // "QUALITY LOCAL REFINEMENT OF TETRAHEDRAL MESHES BASED ON 8 - SUBTETRAHEDRON SUBDIVISION"
    // Figure 2
    std::vector<mesh::LocalIndex> t0_conn = { c[0],
                                              c[4],
                                              c[5],
                                              c[6] }; // reuse current element for t0
    mesh.elements().add(ElementType::T3D4N, { c[4], c[1], c[7], c[8] });
    mesh.elements().add(ElementType::T3D4N, { c[5], c[7], c[2], c[9] });
    mesh.elements().add(ElementType::T3D4N, { c[6], c[8], c[9], c[3] });
    mesh.elements().add(ElementType::T3D4N, { c[4], c[8], c[6], c[5] });
    mesh.elements().add(ElementType::T3D4N, { c[4], c[7], c[8], c[5] });
    mesh.elements().add(ElementType::T3D4N, { c[9], c[5], c[7], c[8] });
    mesh.elements().add(ElementType::T3D4N, { c[9], c[6], c[5], c[8] });
    for (auto i : utility::range(4))
    {
        element.connectivity()[i] = t0_conn[i];
    }
}

void refine_uniform(Mesh& mesh)
{
    bool is_quadratic = false;
    if (!is_type_linear(mesh.elements()[0].type()))
    {
        convert_to_linear(mesh);
        is_quadratic = true;
    }

    std::map<std::size_t, LocalIndex> edge_hash_to_new_node_map;
    auto n_elements_orig = mesh.elements().size();
    for (auto i : utility::range(n_elements_orig))
    {
        if (mesh.elements()[i].type() == ElementType::T3D4N)
        {
            refine_linear_tetrahedra(mesh, mesh.elements()[i], edge_hash_to_new_node_map);
        }
        else
        {
            throw _UTIL_MSG_EXCEPT("Refinement kernel not implemented for element type.  Extend "
                                   "refinement algorithms to support type " +
                                   to_string(mesh.elements()[i].type()));
        }
    }

    if (is_quadratic)
    {
        convert_to_quadratic(mesh);
    }
}
}