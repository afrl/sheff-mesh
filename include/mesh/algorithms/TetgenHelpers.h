#pragma once

#include "mesh/Mesh.h"
#include "second_party_mesh/tetgen/tetgen.h"

namespace mesh
{
/// Creates a mesh that describes the input to tetgen.
Mesh get_mesh_from_tetgen_input(tetgen::tetgenio& in);

/// Creates a mesh that describes the output from tetgen.
Mesh get_mesh_from_tetgen_output(tetgen::tetgenio& out);

/// @brief Creates a tetrahedralization given a surface mesh using TetGen.
/// Currently, this only supports triangular surface meshes, so every element in the surface mesh
/// must be a triangle.
Mesh tetrahedralize_surface(const Mesh& surface_mesh, const bool delaunay = false);

/// @brief Populates an input struct that TetGen can tetrahedralize from a surface mesh.
std::unique_ptr<tetgen::tetgenio> surface_to_tetgen_input(const Mesh& surface_mesh);

/// @brief Calculates the volume of a single tetrahedra (sign indicates inversion)
double tetrahedron_signed_volume(const Eigen::RowVector3d& p1,
                                 const Eigen::RowVector3d& p2,
                                 const Eigen::RowVector3d& p3,
                                 const Eigen::RowVector3d& p4);

struct SegmentedTetrahedralization
{
    Mesh tri_surface;
    Mesh tetrahedralization;
    unsigned int n_regions;
    std::vector<std::vector<mesh::LocalIndex>> region_to_tet_map;
    simple_math::GeneralMatrix region_interior_points;
};

SegmentedTetrahedralization segmented_tetrahedralization(const Mesh& surface_mesh);
}