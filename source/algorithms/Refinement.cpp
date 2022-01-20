#include "mesh/algorithms/Refinement.h"

#include "mesh/algorithms/Manipulation.h"
#include "simple_math/Utility.h"
#include "utility/Hash.h"

using namespace simple_math::point_operators;

namespace mesh
{

struct RefinementEntry
{
    std::vector<LocalIndex> orig_nodes;
    LocalIndex new_node;
};

struct ElementRefinementEntry
{
    LocalIndex orig_element;
    std::vector<LocalIndex> new_elements;
};

inline std::vector<LocalIndex> create_new_nodes(
  const ElementConstAccessor& element,
  const std::vector<std::vector<ConnIndex>>& split_nodes,
  Mesh& mesh,
  std::map<std::size_t, RefinementEntry>& edge_face_hash_to_new_node_map)
{
    // Find or create new midside nodes (4-9 in the figure above)
    std::vector<mesh::LocalIndex> c =
      std::vector<mesh::LocalIndex>(element.connectivity().begin(), element.connectivity().end());
    std::vector<mesh::LocalIndex> split_conn;
    for (auto& n : split_nodes)
    {
        split_conn.resize(n.size());
        for (auto i : utility::range(n.size()))
        {
            split_conn[i] = c[n[i]];
        }
        std::sort(split_conn.begin(), split_conn.end());

        const auto segment_hash = utility::hash(split_conn);
        auto find_iter          = edge_face_hash_to_new_node_map.find(segment_hash);
        if (find_iter == edge_face_hash_to_new_node_map.end())
        { // If it doesn't exist, add a midside node
            Vector3D new_coord(0, 0, 0);
            for (auto i : utility::range(n.size()))
            {
                new_coord += mesh.nodes()[c[n[i]]].coord();
            }
            new_coord = new_coord / (double)n.size();
            mesh.nodes().add(simple_math::to_point(new_coord));

            edge_face_hash_to_new_node_map[segment_hash] = { split_conn, mesh.nodes().size() - 1 };
            c.push_back(mesh.nodes().size() - 1);
        }
        else
        {
            c.push_back(find_iter->second.new_node);
        }
    }
    return c;
}

void refine_linear_tetrahedron(
  Mesh& mesh,
  const ElementAccessor& element,
  std::map<std::size_t, RefinementEntry>& edge_face_hash_to_new_node_map,
  std::vector<ElementRefinementEntry>& element_refinements)
{
    // Based on connectivity: (with new nodes labeled)
    //      3
    //     /|\
    //    / | \
    //  6/  |8 \9
    //  /   |   \
    // 0----1----2
    //   4  6  5
    //
    // View of bottom face:
    //      5
    // 0---------2
    //  \       /
    //   \     /
    //   4\   /7
    //     \ /
    //      1

    static const std::vector<std::vector<ConnIndex>> split_nodes = {
        { 0, 1 }, // 4
        { 0, 2 }, // 5
        { 0, 3 }, // 6
        { 1, 2 }, // 7
        { 1 ,3 }, // 8
        { 2, 3 }, // 9
    };

    std::vector<mesh::LocalIndex> c =
      create_new_nodes(element, split_nodes, mesh, edge_face_hash_to_new_node_map);

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
    mesh.elements().add(ElementType::T3D4N, { c[4], c[6], c[8], c[5] });
    mesh.elements().add(ElementType::T3D4N, { c[4], c[8], c[7], c[5] });
    mesh.elements().add(ElementType::T3D4N, { c[9], c[5], c[7], c[8] });
    mesh.elements().add(ElementType::T3D4N, { c[9], c[6], c[5], c[8] });
    for (auto i : utility::range(4))
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

void refine_quadratic_tetrahedron(
  Mesh& mesh,
  const ElementAccessor& element,
  std::map<std::size_t, RefinementEntry>& edge_face_hash_to_new_node_map,
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

    static const std::vector<std::vector<ConnIndex>> split_nodes = {
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

    std::vector<mesh::LocalIndex> c =
      create_new_nodes(element, split_nodes, mesh, edge_face_hash_to_new_node_map);

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
    for (auto i : utility::range_as<mesh::ConnIndex>(t0_conn.size()))
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

void refine_linear_hexahedron(
  Mesh& mesh,
  const ElementAccessor& element,
  std::map<std::size_t, RefinementEntry>& edge_face_hash_to_new_node_map,
  std::vector<ElementRefinementEntry>& element_refinements)
{
    // Original connectivity:
    //     3---------2
    //    /|        /|
    //   / |       / |
    //  7---------6  |
    //  |  0------|--1
    //  | /       | /
    //  |/        |/
    //  4---------5

    // New nodes
    //   Back face        Mid face      Front face
    //  3----10---2     16---19---15    7----24---6
    //  |         |     |         |     |         |
    //  11   12   9     20   21   18    25   26   23
    //  |         |     |         |     |         |
    //  0----8----1     13---17---14    4----22---5

    static const std::vector<std::vector<ConnIndex>> split_nodes = {
        { 0, 1 },       // 8
        { 1, 2 },       // 9
        { 2, 3 },       // 10
        { 3, 0 },       // 11
        { 0, 1, 2, 3 }, // 12

        { 0, 4 }, // 13
        { 1, 5 }, // 14
        { 2, 6 }, // 15
        { 3, 7 }, // 16

        { 13, 14 },         // 17
        { 14, 15 },         // 18
        { 15, 16 },         // 19
        { 16, 13 },         // 20
        { 13, 14, 15, 16 }, // 21

        { 4, 5 },       // 22
        { 5, 6 },       // 23
        { 6, 7 },       // 24
        { 7, 4 },       // 25
        { 4, 5, 6, 7 }, // 26
    };

    std::vector<mesh::LocalIndex> c =
      create_new_nodes(element, split_nodes, mesh, edge_face_hash_to_new_node_map);

    std::vector<mesh::LocalIndex> replacement_0_conn = { c[0],  c[8],  c[12], c[11],
                                                         c[13], c[17], c[21], c[20] };
    mesh.elements().add(ElementType::H3D8N,
                        { c[8], c[1], c[9], c[12], c[17], c[14], c[18], c[21] });
    mesh.elements().add(ElementType::H3D8N,
                        { c[12], c[9], c[2], c[10], c[21], c[18], c[15], c[19] });
    mesh.elements().add(ElementType::H3D8N,
                        { c[11], c[12], c[10], c[3], c[20], c[21], c[19], c[16] });
    mesh.elements().add(ElementType::H3D8N,
                        { c[13], c[17], c[21], c[20], c[4], c[22], c[26], c[25] });
    mesh.elements().add(ElementType::H3D8N,
                        { c[17], c[14], c[18], c[21], c[22], c[5], c[23], c[26] });
    mesh.elements().add(ElementType::H3D8N,
                        { c[21], c[18], c[15], c[19], c[26], c[23], c[6], c[24] });
    mesh.elements().add(ElementType::H3D8N,
                        { c[20], c[21], c[19], c[16], c[25], c[26], c[24], c[7] });
    for (auto i : utility::range_as<mesh::ConnIndex>(replacement_0_conn.size()))
    {
        element.connectivity()[i] = replacement_0_conn[i];
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


void refine_quadratic_hexahedron(
  Mesh& mesh,
  const ElementAccessor& element,
  std::map<std::size_t, RefinementEntry>& edge_face_hash_to_new_node_map,
  std::vector<ElementRefinementEntry>& element_refinements)
{
    // Original connectivity:
    //      6----5----4
    //     /|        /|
    //    / 7       / 3
    //  18----17---16 |
    //   |  0----1-|--2
    //  19 /       15/
    //   |/        |/
    //  12----13---14
    //
    //   Back face
    // 11--------10
    //  |         |
    //  |         |
    //  |         |
    //  8---------9

    // New faces (note out of order)
    //
    //       Face 1
    //  6--25--5--24--4
    //  |      |      |
    // 26      31     23
    //  |      |      |
    //  7--32--28-30--3
    //  |      |      |
    // 27      29     22
    //  |      |      |
    //  0--20--1--21--2
    //
    //       Face 3
    // 11--43--35-42--10
    //  |      |      |
    // 44      48     41
    //  |      |      |
    // 36--49--37-47--34
    //  |      |      |
    // 45      46     40
    //  |      |      |
    //  8--38--33-39--9
    //
    //       Face 2
    // 56------55-----54
    //  |      |      |
    //  |      |      |
    //  |      |      |
    // 57------58-----53
    //  |      |      |
    //  |      |      |
    //  |      |      |
    // 50------51-----52
    //
    //       Face 5
    // 18--64--17-63--16
    //  |      |      |
    // 65      70     62
    //  |      |      |
    // 19--71--67-69--15
    //  |      |      |
    // 66      68     61
    //  |      |      |
    // 12--59--13-60--14
    //
    //       Face 4
    // 78------77-----76
    //  |      |      |
    //  |      |      |
    //  |      |      |
    // 79------80-----75
    //  |      |      |
    //  |      |      |
    //  |      |      |
    // 72------73-----74

    static const std::vector<std::vector<ConnIndex>> split_nodes = {
        // Face 1
        { 0, 1 },       // 20
        { 1, 2 },       // 21
        { 2, 3 },       // 22
        { 3, 4 },       // 23
        { 4, 5 },       // 24
        { 5, 6 },       // 25
        { 6, 7 },       // 26
        { 7, 0 },       // 27
        { 1, 3, 5, 7 }, // 28
        { 1, 28 },      // 29
        { 3, 28 },      // 30
        { 5, 28 },      // 31
        { 7, 28 },      // 32
        // Face 3
        { 8, 9 },         // 33
        { 9, 10 },        // 34
        { 10, 11 },       // 35
        { 11, 8 },        // 36
        { 8, 9, 10, 11 }, // 37
        { 8, 33 },        // 38
        { 33, 9 },        // 39
        { 9, 34 },        // 40
        { 34, 10 },       // 41
        { 10, 35 },       // 42
        { 35, 11 },       // 43
        { 11, 36 },       // 44
        { 36, 8 },        // 45
        { 33, 37 },       // 46
        { 34, 37 },       // 47
        { 35, 37 },       // 48
        { 36, 37 },       // 49
        // Face 2
        { 0, 8 },   // 50
        { 1, 33 },  // 51
        { 2, 9 },   // 52
        { 3, 34 },  // 53
        { 4, 10 },  // 54
        { 5, 35 },  // 55
        { 6, 11 },  // 56
        { 7, 36 },  // 57
        { 28, 37 }, // 58
        // Face 5
        { 12, 13 },         // 59
        { 13, 14 },         // 60
        { 14, 15 },         // 61
        { 15, 16 },         // 62
        { 16, 17 },         // 63
        { 17, 18 },         // 64
        { 18, 19 },         // 65
        { 19, 12 },         // 66
        { 13, 15, 17, 19 }, // 67
        { 13, 67 },         // 68
        { 15, 67 },         // 69
        { 17, 67 },         // 70
        { 19, 67 },         // 71
        // Face 4
        { 8, 12 },  // 72
        { 33, 13 }, // 73
        { 9, 14 },  // 74
        { 34, 15 }, // 75
        { 10, 16 }, // 76
        { 35, 17 }, // 77
        { 11, 18 }, // 78
        { 36, 19 }, // 79
        { 37, 67 }, // 80
    };

    // Find or create new midside nodes (4-9 in the figure above)
    std::vector<mesh::LocalIndex> c =
      create_new_nodes(element, split_nodes, mesh, edge_face_hash_to_new_node_map);

    /*
    static const std::vector<std::array<ConnIndex, 20>> new_conns = {
        { 0, 20, 1, 29, 28, 32, 7, 27, 50, 51, 58, 57, 8, 38, 33, 46, 37, 49, 36, 45 },
        { 1, 21, 2, 22, 3, 30, 28, 29, 51, 52, 53, 58, 33, 39, 9, 40, 34, 47, 37, 46 },
        { 28, 30, 3, 23, 4, 24, 5, 31, 58, 53, 54, 55, 37, 47, 34, 41, 10, 42, 35, 48 },
        { 7, 32, 28, 31, 5, 25, 6, 26, 57, 58, 55, 56, 36, 49, 37, 48, 35, 43, 11, 44 },
        { 8, 38, 33, 46, 37, 49, 36, 45, 72, 73, 80, 79, 12, 59, 13, 68, 67, 71, 19, 66 },
        { 33, 39, 9, 40, 34, 47, 37, 46, 73, 74, 75, 80, 13, 60, 14, 61, 15, 69, 67, 68 },
        { 37, 47, 34, 41, 10, 42, 35, 48, 80, 75, 76, 77, 67, 69, 15, 62, 16, 63, 17, 70 },
        { 36, 49, 37, 48, 35, 43, 11, 44, 79, 80, 77, 78, 19, 71, 67, 70, 17, 64, 18, 65 }
    };
    */

    std::vector<mesh::LocalIndex> replacement_0_conn = { c[0],  c[20], c[1],  c[29], c[28],
                                                         c[32], c[7],  c[27], c[50], c[51],
                                                         c[58], c[57], c[8],  c[38], c[33],
                                                         c[46], c[37], c[49], c[36], c[45] };
    mesh.elements().add(ElementType::H3D20N,
                        { c[1],  c[21], c[2],  c[22], c[3], c[30], c[28], c[29], c[51], c[52],
                          c[53], c[58], c[33], c[39], c[9], c[40], c[34], c[47], c[37], c[46] });
    mesh.elements().add(ElementType::H3D20N,
                        { c[28], c[30], c[3],  c[23], c[4],  c[24], c[5],  c[31], c[58], c[53],
                          c[54], c[55], c[37], c[47], c[34], c[41], c[10], c[42], c[35], c[48] });
    mesh.elements().add(ElementType::H3D20N,
                        { c[7],  c[32], c[28], c[31], c[5],  c[25], c[6],  c[26], c[57], c[58],
                          c[55], c[56], c[36], c[49], c[37], c[48], c[35], c[43], c[11], c[44] });
    mesh.elements().add(ElementType::H3D20N,
                        { c[8],  c[38], c[33], c[46], c[37], c[49], c[36], c[45], c[72], c[73],
                          c[80], c[79], c[12], c[59], c[13], c[68], c[67], c[71], c[19], c[66] });
    mesh.elements().add(ElementType::H3D20N,
                        { c[33], c[39], c[9],  c[40], c[34], c[47], c[37], c[46], c[73], c[74],
                          c[75], c[80], c[13], c[60], c[14], c[61], c[15], c[69], c[67], c[68] });
    mesh.elements().add(ElementType::H3D20N,
                        { c[37], c[47], c[34], c[41], c[10], c[42], c[35], c[48], c[80], c[75],
                          c[76], c[77], c[67], c[69], c[15], c[62], c[16], c[63], c[17], c[70] });
    mesh.elements().add(ElementType::H3D20N,
                        { c[36], c[49], c[37], c[48], c[35], c[43], c[11], c[44], c[79], c[80],
                          c[77], c[78], c[19], c[71], c[67], c[70], c[17], c[64], c[18], c[65] });
    for (auto i : utility::range_as<mesh::ConnIndex>(replacement_0_conn.size()))
    {
        element.connectivity()[i] = replacement_0_conn[i];
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
    std::map<std::size_t, RefinementEntry> edge_face_hash_to_new_node_map;
    auto n_elements_orig = mesh.elements().size();
    for (auto i : utility::range(n_elements_orig))
    {
        if (mesh.elements()[i].type() == ElementType::T3D4N)
        {
            refine_linear_tetrahedron(mesh,
                                      mesh.elements()[i],
                                      edge_face_hash_to_new_node_map,
                                      element_refinements);
        }
        else if (mesh.elements()[i].type() == ElementType::T3D10N)
        {
            refine_quadratic_tetrahedron(mesh,
                                         mesh.elements()[i],
                                         edge_face_hash_to_new_node_map,
                                         element_refinements);
        }
        else if (mesh.elements()[i].type() == ElementType::H3D8N)
        {
            refine_linear_hexahedron(mesh,
                                     mesh.elements()[i],
                                     edge_face_hash_to_new_node_map,
                                     element_refinements);
        }
        else if (mesh.elements()[i].type() == ElementType::H3D20N)
        {
            refine_quadratic_hexahedron(mesh,
                                        mesh.elements()[i],
                                        edge_face_hash_to_new_node_map,
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
        for (auto& refined_edge : edge_face_hash_to_new_node_map)
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