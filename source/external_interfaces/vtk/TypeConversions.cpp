#include "mesh/external_interfaces/vtk/TypeConversions.h"

#include <fstream>
#include <string>

namespace mesh
{
int get_vtk_element_type(const ElementType& geom_type)
{
    // 1D
    if (geom_type == ElementType::L1D2N)
        return 3;
    else if (geom_type == ElementType::L1D3N)
        return 21;
    // 2D
    else if (geom_type == ElementType::T2D3N)
        return 5;
    else if (geom_type == ElementType::Q2D4N)
        return 9;
    else if (geom_type == ElementType::T2D6N)
        return 22;
    else if (geom_type == ElementType::Q2D8N)
        return 23;
    else if (geom_type == ElementType::POLYGON)
        return 7;
    // else if (geom_type == ElementType::POLYGON)
    //	return 7;
    // 3D
    else if (geom_type == ElementType::T3D4N)
        return 10;
    else if (geom_type == ElementType::P3D5N)
        return 14;
    else if (geom_type == ElementType::W3D6N)
        return 13;
    else if (geom_type == ElementType::H3D8N)
        return 12;
    else if (geom_type == ElementType::T3D10N)
        return 24;
    else if (geom_type == ElementType::P3D13N)
        return 27;
    else if (geom_type == ElementType::W3D15N)
        return 26;
    else if (geom_type == ElementType::H3D20N)
        return 25;
    else if (geom_type == ElementType::POLYHEDRON)
        return 42;
    else
        throw _UTIL_MSG_EXCEPT("VTK type conversion not implemented.");
}

std::vector<ConnIndex> get_vtk_permutation(const ElementType& geom_type, const ConnIndex& n_nodes)
{
    std::vector<ConnIndex> permutation(n_nodes);
    // Fill permutation with 0 - num_nodes-1
    std::iota(permutation.begin(), permutation.end(), 0);

    // 1D
    switch (geom_type)
    {
        case ElementType::L1D2N: return permutation; // Unchanged
        case ElementType::L1D3N:
            throw _UTIL_MSG_EXCEPT("L1D2N type not implemented.");
            // 2D
        case ElementType::T2D3N: return permutation;                // Unchanged
        case ElementType::Q2D4N: return permutation;                // Unchanged
        case ElementType::T2D6N: return { 0, 2, 4, 1, 3, 5 };       // Different
        case ElementType::Q2D8N: return { 0, 2, 4, 6, 1, 3, 5, 7 }; // Different
        case ElementType::POLYGON:
            return permutation; // Unchanged
        // 3D
        case ElementType::T3D4N: return permutation;          // Unchanged
        case ElementType::P3D5N: return permutation;          // Unchanged
        case ElementType::W3D6N: return { 0, 2, 1, 3, 5, 4 }; // Different
        case ElementType::H3D8N: return permutation;          // Unchanged
        case ElementType::T3D10N: return permutation;         // Unchanged
        case ElementType::P3D13N: return permutation;         // Unchanged
        case ElementType::W3D15N:
            return { 0, 4, 2, 9, 13, 11, 5, 3, 1, 14, 12, 10, 6, 8, 7 }; // Different
        case ElementType::H3D20N:
            return {
                0, 2, 4, 6, 12, 14, 16, 18, 1, 3, 5, 7, 13, 15, 17, 19, 8, 9, 10, 11
            };                                            // Different
        case ElementType::POLYHEDRON: return permutation; // Unchanged
        default: return {};
    }
}
}