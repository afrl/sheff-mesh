#include "mesh/algorithms/Refinement.h"

#include "mesh/algorithms/Manipulation.h"
#include "simple_math/Utility.h"
#include "utility/Hash.h"

using namespace simple_math::point_operators;

namespace mesh
{

struct EdgeRefinementEntry
{
    std::array<LocalIndex, 2> orig_nodes;
    LocalIndex new_node;
};

struct ElementRefinementEntry
{
    LocalIndex orig_element;
    std::vector<LocalIndex> new_elements;
};

void refine_linear_tetrahedra(Mesh& mesh,
                              const ElementAccessor& element,
                              std::map<std::size_t, EdgeRefinementEntry>& edge_hash_to_new_node_map,
                              std::vector<ElementRefinementEntry>& element_refinements)
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
            edge_hash_to_new_node_map[edge_hash] = { { edge_conn[0], edge_conn[1] },
                                                     mesh.nodes().size() - 1 };
            c.push_back(mesh.nodes().size() - 1);
        }
        else
        {
            c.push_back(find_iter->second.new_node);
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

    std::vector<LocalIndex> new_elements;
    for (auto i : utility::range(mesh.elements().size() - 8, mesh.elements().size()))
    {
        new_elements.push_back(mesh.elements().size() - 1);
    }
    element_refinements.push_back({ element.local_index(), new_elements });
}

void refine_quadratic_tetrahedra(
  Mesh& mesh,
  const ElementAccessor& element,
  std::map<std::size_t, EdgeRefinementEntry>& edge_hash_to_new_node_map,
  std::vector<ElementRefinementEntry>& element_refinements)
{
    // Original connectivity:
    //         3
    //        /|\
    //       / | \
    //      /  |  \
    //    7/   |8  \9
    //    /    |    \
    //   /     |     \
    //  /      |      \
    // 0---4---1---5---2
    //
    // View of bottom face:
    //  0-----6-----2
    //   \         /
    //    \       /
    //    4\     /5
    //      \   /
    //       \ /
    //        1

    static const std::vector<std::array<ConnIndex, 2>> segments_to_split = {
        { 0, 4 }, // 10
        { 4, 1 }, // 11
        { 1, 5 }, // 12
        { 5, 2 }, // 13
        { 2, 6 }, // 14
        { 6, 0 }, // 15
        { 0, 7 }, // 16
        { 7, 3 }, // 17
        { 1, 8 }, // 18
        { 8, 3 }, // 19
        { 2, 9 }, // 20
        { 9, 3 }, // 21
        { 4, 5 }, // 22
        { 5, 6 }, // 23
        { 6, 4 }, // 24
        { 4, 8 }, // 25
        { 8, 7 }, // 26
        { 7, 4 }, // 27
        { 5, 9 }, // 28
        { 9, 8 }, // 29
        { 8, 5 }, // 30
        { 6, 7 }, // 31
        { 7, 9 }, // 32
        { 9, 6 }, // 33
        { 8, 6 }  // 34
    };

    auto orig_conn = element.connectivity();
    auto edge_info = get_edge_info(element.type());

    // Find or create new midside nodes (4-9 in the figure above)
    std::vector<mesh::LocalIndex> c =
      std::vector<mesh::LocalIndex>(orig_conn.begin(), orig_conn.end());
    for (auto& segment : segments_to_split)
    {
        const auto segment_hash =
          utility::hash_sorted(orig_conn[segment[0]], orig_conn[segment[1]]);
        auto find_iter = edge_hash_to_new_node_map.find(segment_hash);
        if (find_iter == edge_hash_to_new_node_map.end())
        { // If it doesn't exist, add a midside node
            auto new_coord = (mesh.nodes()[orig_conn[segment[0]]].coord() +
                              mesh.nodes()[orig_conn[segment[1]]].coord()) /
                             2.0;
            mesh.nodes().add(simple_math::to_point(new_coord));
            edge_hash_to_new_node_map[segment_hash] = { { orig_conn[segment[0]],
                                                          orig_conn[segment[1]] },
                                                        mesh.nodes().size() - 1 };
            c.push_back(mesh.nodes().size() - 1);
        }
        else
        {
            c.push_back(find_iter->second.new_node);
        }
    }

    // Create the new tetrahedra
    // See https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.308.2549&rep=rep1&type=pdf
    // "QUALITY LOCAL REFINEMENT OF TETRAHEDRAL MESHES BASED ON 8 - SUBTETRAHEDRON SUBDIVISION"
    // Figure 2
    // reuse current element for t0
    std::vector<mesh::LocalIndex> t0_conn = { c[0],  c[4],  c[6],  c[7],  c[10],
                                              c[24], c[15], c[16], c[27], c[31] };
    mesh.elements().add(ElementType::T3D10N,
                        { c[4], c[1], c[5], c[8], c[11], c[12], c[22], c[25], c[18], c[30] });
    mesh.elements().add(ElementType::T3D10N,
                        { c[5], c[2], c[6], c[9], c[13], c[14], c[23], c[28], c[20], c[33] });
    mesh.elements().add(ElementType::T3D10N,
                        { c[7], c[8], c[9], c[3], c[26], c[29], c[32], c[17], c[19], c[21] });
    mesh.elements().add(ElementType::T3D10N,
                        { c[4], c[8], c[6], c[7], c[25], c[34], c[24], c[27], c[26], c[31] });
    mesh.elements().add(ElementType::T3D10N,
                        { c[4], c[8], c[5], c[6], c[25], c[30], c[22], c[24], c[34], c[23] });
    mesh.elements().add(ElementType::T3D10N,
                        { c[6], c[5], c[9], c[8], c[23], c[28], c[33], c[34], c[30], c[29] });
    mesh.elements().add(ElementType::T3D10N,
                        { c[9], c[7], c[6], c[8], c[32], c[31], c[33], c[29], c[26], c[34] });
    for (auto i : utility::range(t0_conn.size()))
    {
        element.connectivity()[i] = t0_conn[i];
    }

    std::vector<LocalIndex> new_elements;
    for (auto i : utility::range(mesh.elements().size() - 7, mesh.elements().size()))
    {
        mesh.elements()[i].local_coord_sys().copy(element.local_coord_sys());
        mesh.elements()[i].component_id() = element.component_id();
        new_elements.push_back(i);
    }
    element_refinements.push_back({ element.local_index(), new_elements });
}

void refine_uniform(Mesh& mesh)
{
    std::vector<ElementRefinementEntry> element_refinements;
    std::map<std::size_t, EdgeRefinementEntry> edge_hash_to_new_node_map;
    auto n_elements_orig = mesh.elements().size();
    for (auto i : utility::range(n_elements_orig))
    {
        if (mesh.elements()[i].type() == ElementType::T3D4N)
        {
            refine_linear_tetrahedra(mesh,
                                     mesh.elements()[i],
                                     edge_hash_to_new_node_map,
                                     element_refinements);
        }
        else if (mesh.elements()[i].type() == ElementType::T3D10N)
        {
            refine_quadratic_tetrahedra(mesh,
                                        mesh.elements()[i],
                                        edge_hash_to_new_node_map,
                                        element_refinements);
        }
        else
        {
            throw _UTIL_MSG_EXCEPT("Refinement kernel not implemented for element type.  Extend "
                                   "refinement algorithms to support type " +
                                   to_string(mesh.elements()[i].type()));
        }
    }

    // Carry over the node groups
    std::vector<bool> is_node_in_group(mesh.nodes().size());
    for (auto& group : mesh.nodes().groups())
    {
        std::fill(is_node_in_group.begin(), is_node_in_group.end(), false);
        for (auto i : group.node_indices)
        {
            is_node_in_group[i] = true;
        }
        for (auto& refined_edge : edge_hash_to_new_node_map)
        {
            if (is_node_in_group[refined_edge.second.orig_nodes[0]] &&
                is_node_in_group[refined_edge.second.orig_nodes[1]])
            {
                group.node_indices.push_back(refined_edge.second.new_node);
            }
        }
    }

    // Carry over the element groups
    std::vector<bool> is_element_in_group(mesh.elements().size());
    for (auto& group : mesh.elements().groups())
    {
        std::fill(is_element_in_group.begin(), is_element_in_group.end(), false);
        for (auto i : group.element_indices)
        {
            is_element_in_group[i] = true;
        }
        for (auto& refined_element : element_refinements)
        {
            if (is_element_in_group[refined_element.orig_element])
            {
                for (auto& e_index : refined_element.new_elements)
                {
                    if (e_index != refined_element.orig_element)
                    {
                        group.element_indices.push_back(e_index);
                    }
                }
            }
        }
    }
}
}