#pragma once

#include "mesh/Mesh.h"
#include "mesh/algorithms/Search.h"
#include "mesh/kd_tree/KDTree.h"

#include <utility>

namespace mesh
{
/// @brief A simple struct that holds the result of find_surface_quadrature()
struct SurfaceQuadratureInfo
{
    /// Index of the volume elements that contains each surface quadrature point.
    /// Size: # of quadrature points
    std::vector<LocalIndex> parent_volume_elements;
    /// Coordinate of each surface quadrature point with respect to the parametric coordinate system
    /// of the volume elements that contains point (see parent_volume_elements).
    std::vector<Point> parent_element_parametric_coordinates;
    /// Coordinate of each surface quadraturepoint in the global coordinate system.
    std::vector<Point> global_coordinates;
    /// Outward normal vector for each surface quadrature point expressed in the global coordinate
    /// system.
    std::vector<simple_math::Vector3D> normals;
    /// Jacobian (of surface element
    std::vector<double> jxw_terms;
};

/// @brief Finds the Gauss quadrature points and other relevant information for all surface elements
/// whose nodes are a subset of the given node indices.
///
/// @param extra_quadrature_order - the quadrature is automatically selected to integrate
///        a constant function over the face.   This indicates how many extra orders or add.
///        For a linear 4-node quad, only 1 quad point is needed, so if extra_quadrature_order
///        is 1, then you will end up with a 2nd order quadrature (2x2).
SurfaceQuadratureInfo find_surface_quadrature(
  const Mesh& mesh,
  const BoundaryInfo& boundary_info,
  const std::vector<LocalIndex>& surface_node_indices,
  const simple_math::QuadratureOrder& extra_quadrature_order = 0);

/// @brief A simple struct to hold the result of a call to find_connections
struct ConnectionResult
{
    /// Number of of elements each point is connected to.  If a point lies on the boundary of
    /// two or more elements, then it can be more than 1.  If the point is too far away from any
    /// element, then it will be 0.
    std::vector<ConnIndex> n_connected_elements_by_point;
    /// An unrolled list of element indices to which the points are connected.  For example, if
    /// there are 3 points with the following connections:
    /// 0 -> { element 0, element 2 }
    /// 1 -> none
    /// 2 -> { element 1 }
    /// Then the vector will contain { element 0, element 2, element 3 }
    /// Note: n_connected_elements_by_point would contain { 2, 0, 1 }
    std::vector<LocalIndex> connected_elements;
    /// Following the order of connected_elements, this holds the parametric coordiante of the
    /// corresponding given point in parametric coordinate system of the respective connected volume
    /// element.
    std::vector<Point> parametric_coordinates_in_connected_elements;
};

/// @brief Finds any connections between a list of points and any elements in a mesh.
///
/// @param search_radius - the radius around each point that an element will be considered (this
///        should be the largest edge in the mesh ideally)
/// @param match_tolerance - if a point does not lie in any element, then this tolerance is used to
/// determine if the point is close enough to be connected to the respective element.  The tolerance
/// refers to the parametric coordinate of the element.  In other words, the point is "close enough"
/// if each component its parametric coordinate is within (-1 - match_tolerance) and
/// (1 + match_tolerance).  Also, the tolerance is not applicable when a point is within any
/// element (even if it would be within tolerance to a different element).
///
/// Different situations that can occur:
/// element 0 element 1
/// o--------o--------o
/// |        |        |
/// |        x p1     |
/// |        |        |
/// o---x----o--------o
///    p0        x p2           x p3
/// p0 is connected to element 0, since it lies on it alone.
/// p1 is connected to both elements 0 and 1, since it lies on the face of both
/// p2 is connected to element 1, since it parametric coordinate is very close to 1 (assuming a
///    lenient tolerance)
/// p3 is not connected to any element.
ConnectionResult find_connections(const Mesh& mesh,
                                  const NodeContainerKDTree& element_centroids_kd_tree,
                                  const std::vector<simple_math::Point>& points_to_connect,
                                  const double search_radius,
                                  const double match_tolerance);

// Do we need this?
/*
void find_connections(
  const Mesh& mesh,
  const std::vector<std::unique_ptr<mesh::NodeContainerKDTree>>&
    element_centroids_kd_tree_by_component,
  const std::vector<std::array<simple_math::Point, 2>>& bounding_box_by_component,
  const std::vector<ComponentIndex>& components_to_connect_to,
  const std::vector<simple_math::Point>& points_to_connect);
*/
}