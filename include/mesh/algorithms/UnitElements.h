#pragma once

#include "mesh/Mesh.h"

namespace mesh
{
/// Given a reference to a mesh, adds a linear hex element over the extent of -1 to 1 (3D)
static Mesh get_unit_linear_hex()
{
    Mesh mesh;
    mesh.nodes().add({ -1, -1, -1 });
    mesh.nodes().add({ -1, 1, -1 });
    mesh.nodes().add({ -1, 1, 1 });
    mesh.nodes().add({ -1, -1, 1 });
    mesh.nodes().add({ 1, -1, -1 });
    mesh.nodes().add({ 1, 1, -1 });
    mesh.nodes().add({ 1, 1, 1 });
    mesh.nodes().add({ 1, -1, 1 });

    mesh.elements().add(mesh::ElementType::H3D8N, { 0, 1, 2, 3, 4, 5, 6, 7 });

    return mesh;
}

static Mesh get_unit_P3D5N()
{
    Mesh mesh;
    mesh.nodes().add({ -1, -1, -1 });
    mesh.nodes().add({ 1, -1, -1 });
    mesh.nodes().add({ 1, 1, -1 });
    mesh.nodes().add({ -1, 1, -1 });
    mesh.nodes().add({ 0, 0, 1 });

    mesh.elements().add(mesh::ElementType::P3D5N, { 0, 1, 2, 3, 4 });

    return mesh;
}
}