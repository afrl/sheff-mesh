#include "mesh/kd_tree/KDTree.h"

using namespace mesh;

int main()
{
    NodeContainer nodes;
    nodes.add({ 0, 0, 0 });
    nodes.add({ 1, 1, 1 });

    NodeContainerKDTree tree(nodes);
    auto n = tree.closest_node({ 0.75, 0.5, 0.5 });
    assert(n.local_index() == 1);

    return 0;
}