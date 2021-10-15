#include "mesh/algorithms/Connections.h"

#include "mesh/algorithms/IO.h"
#include "mesh/algorithms/IsoparametricMappingUtils.h"
#include "mesh/basis/BasisFunctionsCache.h"
#include "mesh/common/QuadratureHelpers.h"
#include "simple_math/Utility.h"
#include "simple_math/quadrature/Quadrature.h"

#include <unordered_set>

namespace mesh
{
SurfaceQuadratureInfo find_surface_quadrature(
  const Mesh& mesh,
  const BoundaryInfo& boundary_info,
  const std::vector<LocalIndex>& surface_node_indices,
  const simple_math::QuadratureOrder& extra_quadrature_order)
{
    // Build a set to identify the surface elements more quickly
    std::set<LocalIndex> surface_nodes;
    for (auto& i : surface_node_indices)
    {
        surface_nodes.insert(i);
    }

    // Find the subset of boundary elements that compose the surface
    // that will connected to the secondary cluster.
    std::vector<mesh::LocalIndex> surface_elements;
    for (auto e : boundary_info.surface)
    {
        bool all_nodes_on_surfaces = true;
        for (auto n_id : e.connectivity())
        {
            if (surface_nodes.find(n_id) == surface_nodes.end())
            {
                all_nodes_on_surfaces = false;
                break;
            }
        }
        if (all_nodes_on_surfaces)
        {
            surface_elements.push_back(e.local_index());
        }
    }

    // mesh::write_vtk(boundary_info.surface.mesh(), "surface.vtk", false);

    SurfaceQuadratureInfo result;
    // Simply a good guess at how many might be needed (typically needs one reallocation at most).
    result.global_coordinates.reserve(4 * surface_elements.size());
    result.jxw_terms.reserve(4 * surface_elements.size());
    result.normals.reserve(4 * surface_elements.size());
    result.parent_element_parametric_coordinates.reserve(4 * surface_elements.size());
    result.parent_volume_elements.reserve(4 * surface_elements.size());

    BasisFunctionsCache basis_cache;
    for (auto s : utility::range(surface_elements.size()))
    {
        auto& surface_e_id = surface_elements[s];
        auto surface_e     = boundary_info.surface[surface_e_id];

        auto surface_e_quadrature_scheme = select_gauss_quadrature_scheme(surface_e.type());
        surface_e_quadrature_scheme.integration_order += extra_quadrature_order;
        simple_math::Quadrature quadrature(surface_e_quadrature_scheme, nullptr);
        basis_cache.setup(surface_e.type(), quadrature);
        auto& quad_to_basis_values_map      = basis_cache.values();
        auto& quad_to_basis_derivatives_map = basis_cache.derivatives();

        auto quad_global_coords =
          quadrature_parametric_to_global(surface_e, quadrature, basis_cache);
        for (auto i : utility::range(quad_global_coords.rows()))
        {
            result.global_coordinates.push_back(simple_math::to_point(quad_global_coords.row(i)));
        }

        auto normal = calc_outward_normal(surface_e);
        // Assume the same normal for every quadrature point
        for (auto i : utility::range(quadrature.n_quad_points()))
        {
            result.normals.push_back(normal);
        }

        auto surface_e_nodal_coords = get_nodal_global_coords(surface_e);
        for (auto q : utility::range(quadrature.n_quad_points()))
        {
            // This is a 2x3 matrix that represents:
            // \left| \frac{ \partial \bar{g}_i^\beta (\bar{\alpha}) }{ \bar{\alpha } \right|
            // where \bar{g}_i^\beta (\bar{\alpha}) is the map from the parametric coordinate system
            // of the triangle, \bar{\alpha}, to the parametric coordinate system of the parent
            // element.
            Eigen::Matrix<double, 2, 3, Eigen::RowMajor> tri_to_element_jacobian =
              quad_to_basis_derivatives_map[q].transpose() * surface_e_nodal_coords;
            auto tri_to_element_jacobian_det =
              tri_to_element_jacobian.row(0).cross(tri_to_element_jacobian.row(1)).norm();

            result.jxw_terms.push_back(tri_to_element_jacobian_det * quadrature.weight(q));
        }

        for (auto q : utility::range(quadrature.n_quad_points()))
        {
            auto volume_e_id = boundary_info.parent_volume_elements[surface_e_id];
            result.parent_volume_elements.push_back(volume_e_id);
            auto volume_e = mesh.elements()[volume_e_id];

            auto param_result = global_to_parametric(volume_e, quadrature.parametric_point(q));
            // Since this is finding the parametric coordinate of a point on the surface of the
            // element, the search should always be succeessful
            assert(param_result.first);

            result.parent_element_parametric_coordinates.push_back(
              simple_math::to_point(param_result.second));
        }
    }

    return result;
}

ConnectionResult find_connections(const Mesh& mesh,
                                  const NodeContainerKDTree& element_centroids_kd_tree,
                                  const std::vector<simple_math::Point>& points_to_connect,
                                  const double search_radius,
                                  const double search_tolerance)
{
    // mesh::write_vtk(mesh, "mesh.vtk", false);

    ConnectionResult result;
    result.n_connected_elements_by_point.resize(points_to_connect.size(), 0);
    result.connected_elements.reserve(points_to_connect.size());
    result.parametric_coordinates_in_connected_elements.reserve(points_to_connect.size());

    for (auto i : utility::indices(points_to_connect))
    {
        auto& p = points_to_connect[i];

        // Find elements that the point lies in or the closest element within the tolerance
        auto closest_element_centroids = element_centroids_kd_tree.in_radius(p, search_radius);
        bool found_element             = false;
        LocalIndex closest_element_id  = 0;
        double closest_distance        = std::numeric_limits<double>::max();
        for (auto e_centroid : closest_element_centroids)
        {
            auto e               = mesh.elements()[e_centroid.local_index()];
            auto contains_result = contains_point2(e, p);

            if (contains_result.first) // contains the point
            {
                result.n_connected_elements_by_point[i] += 1;
                result.connected_elements.push_back(e_centroid.local_index());

                auto param_result = global_to_parametric(e, p);
                // Since this is finding the parametric coordinate of a point in the element,
                // the search should always be succeessful
                assert(param_result.first);
                result.parametric_coordinates_in_connected_elements.push_back(
                  simple_math::to_point(param_result.second));

                found_element = true;
            }

            if (contains_result.second < closest_distance) // track closest in case none contain it
            {
                closest_element_id = e_centroid.local_index();
                closest_distance   = contains_result.second;
            }
        }

        if (!found_element && closest_distance != std::numeric_limits<double>::max())
        {
            auto closest_e    = mesh.elements()[closest_element_id];
            auto param_result = mesh::global_to_parametric(closest_e, p);
            if (param_result.first && param_result.second[0] > -1.0 - search_tolerance &&
                param_result.second[1] > -1.0 - search_tolerance &&
                param_result.second[2] > -1.0 - search_tolerance &&
                param_result.second[0] < 1.0 + search_tolerance &&
                param_result.second[1] < 1.0 + search_tolerance &&
                param_result.second[2] < 1.0 + search_tolerance)
            {
                result.n_connected_elements_by_point[i] = 1;
                result.connected_elements.push_back(closest_element_id);
                result.parametric_coordinates_in_connected_elements.push_back(
                  simple_math::to_point(param_result.second));
            }
        }
    }

    return result;
}
}