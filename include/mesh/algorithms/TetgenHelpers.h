#pragma once

#include "mesh/Mesh.h"
#include "second_party_mesh/tetgen/tetgen.h"

namespace mesh
{
/// Extracts surface mesh from tetgenio.
Mesh get_surface_from_tetgenio(tetgen::tetgenio& io);

/// Extracts volume mesh from tetgenio.
Mesh get_volume_from_tetgenio(tetgen::tetgenio& io);

/// @brief Populates an input struct that TetGen can tetrahedralize from a surface mesh.
std::unique_ptr<tetgen::tetgenio> surface_to_tetgen_input(const Mesh& surface_mesh);

/// @brief Creates a tetrahedralization given a surface mesh using TetGen.
/// Currently, this only supports triangular surface meshes, so every element in the surface mesh
/// must be a triangle.
Mesh tetrahedralize_surface(const Mesh& surface_mesh, const bool delaunay = false);

/// @brief Populates an input struct that TetGen can tetrahedralize from a volume mesh.
std::unique_ptr<tetgen::tetgenio> volume_to_tetgen_input(const Mesh& volume_mesh);

/// @brief Parameters for algorithms that improve tetrahedralizations
struct ImproveTetParameters
{
    /// Target ratio of a tetrahedron edge length to the diagonal of the bounding box of the given
    /// mesh. A value of -1 indicates there is no size constraint. Smaller positive values will
    /// result in a denser mesh. For reference, 1e-1 will result in a very coarse mesh, while 1e-2
    /// will result in a relatively fine mesh.  Values smaller than 1e-2 can result in very large
    /// meshes.
    ///
    /// Note: for TetGen this is used to calculate the maximum element volume, which is the quality
    /// control for that algorithm.
    double rel_target_edge_length = -1.0;
    /// Maximum value of the ratio of the radius of a sphere circumscribing a tetrahedron to the
    /// minimum edge length of the tetrahedron. Suggested values would be in the range of 1.1 to 1.5
    /// for very high quality. A value closer to 1 will result in higher quality but also more
    /// tetrahedra (most likely). A value of 0 will indicate that there is no maximum.
    double max_radius_edge_ratio = 1.2;
    /// Minimum dihedral angle between any two faces of a tetrahedron. Specifically, there are 6
    /// dihedral angles for a tetrahedron. For a perfectly "regular" tetrahedron, each angle would
    /// be 60 degrees. Realistically, a value of 60 would never converge. Suggested values are in
    /// the range of 3.0 to 10.0.
    double min_dihedral_angle = 5.0;
    /// Maximum number of points that can be added during mesh improvement. A value of -1 indicates
    /// there is no limit. A value of 0 results in no change to the mesh. A value > 1 will result in
    /// bounding the resulting # of vertices at the cost of not losing quality guarantees.
    int max_inserted_points = -1;
};

/// @brief Improves a given a tetrahedralization using TetGen.
///
/// Note: currently only supports 4-node tetrahedra as an input element type.
/// @see ImproveTetParameters
Mesh improve_tetrahedralization(const Mesh& volume_mesh, const ImproveTetParameters& params);

/// @brief Improves a given a tetrahedralization using TetGen.
///
/// Same as other function with same name, but requires less memory if the tetgetio struct is handy.
/// Note: currently only supports 4-node tetrahedra as an input element type.
/// @see ImproveTetParameters
Mesh improve_tetrahedralization(tetgen::tetgenio& input, const ImproveTetParameters& params);

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