#include "mesh/algorithms/Creation.h"
#include "mesh/algorithms/IO.h"
#include "mesh/relations/ElementAdjacencyList.h"

using namespace mesh;

int main()
{
    Mesh linear_grid = make_linear_grid({ 0, 1, 2 }, { 0, 1, 2 }, { 0, 1, 2 });

    ElementAdjacencyList adj;
    adj.update(linear_grid);

    for (auto i : utility::range(linear_grid.elements().size()))
    {
        assert(adj.get_adjacent_element_indices(i, ElementAdjacencyType::Face).size() == 3);
        assert(adj.get_adjacent_element_indices(i, ElementAdjacencyType::FaceEdge).size() == 6);
        assert(adj.get_adjacent_element_indices(i, ElementAdjacencyType::FaceEdgeNode).size() == 7);
    }

    Mesh mixed_grid = make_mixed_element_grid({ 0, 0, 0 }, { 1, 1, 1 }, { 1, 1, 1 });
    adj.update(mixed_grid);
    std::vector<LocalIndex> adjacency_markers(mixed_grid.elements().size());
    for (auto e : mixed_grid.elements())
    {
        for (auto& v : adjacency_markers)
        {
            v = 0;
        }

        adjacency_markers[e.local_index()] = 1;

        for (auto n : adj.get_adjacent_element_indices(e.local_index(), ElementAdjacencyType::Face))
        {
            adjacency_markers[n] = 2;
        }

        write_vtk(mixed_grid,
                  "mixed_adj_" + std::to_string(e.local_index()) + ".vtk",
                  false,
                  {},
                  {},
                  {},
                  adjacency_markers);
    }

    return 0;
}