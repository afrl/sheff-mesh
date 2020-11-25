#pragma once

#include "mesh/Mesh.h"

namespace mesh
{
/// @brief A struct to store information about a tow
struct VTMSTowMesh
{
    enum class MeshType
    {
        POLY,
        TRI,
        TET
    };

    enum class NormalDirection
    {
        INWARD,
        OUTWARD
    };

    enum class SurfaceElementType
    {
        INTERIOR, ///< Polygon or triangle that is interior to the textile model (not clipped)
        BOUNDARY  ///< Polygon or triangle that lies on the boundary of the textile model (clipped)
    };

    int material_number;
    std::vector<Point> medials;
    Mesh mesh;
    std::vector<SurfaceElementType> surface_element_types;
    NormalDirection surface_normal_direction;
    MeshType type;
};

/// @brief Reads a .ctw file from VTMS into a mesh as a collection of polygonal elements that define
/// the bounds of a tow.
/// @param file_path - relative (to cwd) or absolute path to .ctw file
VTMSTowMesh read_vtms_ctw(const std::string& file_path);

/// @brief Takes a polygonal discretization of a tow's surface and transforms it into a
/// quality triangulation of the tow's surface.
/// @param tow_surface - polygonal surface mesh for a tow returned by read_vtms_ctw
void triangulate_tow_surface(VTMSTowMesh& tow_surface);

/// @brief Write a triangulated tow surface to a .ply file, which can be read in by TetWild.
void write_to_ply(const Mesh& triangular_surface, const std::string& file_path);

Mesh create_triangular_matrix_box(const Point& neg_bound, const Point& pos_bound);

/// @brief Takes a triangular discretization of a tow's surface and transforms it into a
/// quality tetrahedralization.
/// @param tow_triangular_surface - triangular surface mesh for a tow returned by
/// triangulate_tow_surface
void tetrahedralize_tow(const VTMSTowMesh& tow_triangular_surface);
}