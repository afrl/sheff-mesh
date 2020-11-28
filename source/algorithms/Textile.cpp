#include "mesh/algorithms/Textile.h"

#include "mesh/algorithms/IO.h"
#include "mesh/algorithms/Search.h"
#include "mesh/relations/NodeToElementMap.h"
#include "second_party_mesh/triangle/triangle.h"
#include "simple_math/Utility.h"
#include "utility/ASCIIStreamHelpers.h"
#include "utility/Hash.h"

#include <fstream>
#include <iomanip>
#include <type_traits>

using namespace simple_math::point_operators;

namespace mesh
{
//////////////////////////////////////////////////////////////////////////////////////////////////
// Private helpers

/// Determines which coordinate axis aligns most closely with the tow path.
/// @returns 0-based index of axis
int8_t find_axis_closest_aligned_along_tow(const VTMSTowMesh& tow_mesh)
{
    auto diff = tow_mesh.medials.front() - tow_mesh.medials.back();
    if (diff[0] > diff[1] && diff[0] > diff[2])
    {
        return 0;
    }
    else if (diff[1] > diff[2])
    {
        return 1;
    }
    else
    {
        return 2;
    }
}

std::vector<std::pair<ConnIndex, ConnIndex>> polygon_segments(const ElementConstAccessor& element)
{
    std::vector<std::pair<ConnIndex, ConnIndex>> segments;
    auto& conn = element.connectivity();
    for (auto i = 0; i < conn.size() - 1; i++)
    {
        segments.push_back({ i, i + 1 });
    }
    segments.push_back({ static_cast<ConnIndex>(conn.size() - 1), 0 });
    return segments;
}

// Implements Newell's method to calculate a polygons normal vector
simple_math::Vector3D calculate_polygon_unit_normal(const ElementConstAccessor& polygon_element)
{
    simple_math::Vector3D normal = simple_math::Vector3D::Zero();
    for (auto segment : polygon_segments(polygon_element))
    {
        auto p0 = polygon_element.node(segment.first).coord();
        auto p1 = polygon_element.node(segment.second).coord();
        normal += p0.cross(p1);
    }
    return normal.normalized();
}

void snap_boundary_points(VTMSTowMesh& tow_surface)
{
    _UTIL_ASSERT(tow_surface.surface_element_types.size() == tow_surface.mesh.elements().size(),
                 "The surface element types must match the number of elements in the mesh.");

    auto extent = find_extent(tow_surface.mesh);

    for (auto element_index :
         utility::range_as<LocalIndex>(tow_surface.surface_element_types.size()))
    {
        auto& element = tow_surface.mesh.elements()[element_index];
        auto& type    = tow_surface.surface_element_types[element_index];

        if (type == VTMSTowMesh::SurfaceElementType::BOUNDARY)
        {
            // Find which plane we are moving this to
            auto normal     = calculate_polygon_unit_normal(element);
            int normal_axis = -1;
            for (auto i = 0; i < normal.size(); i++)
            {
                if (fabs(normal[i]) > 0.9)
                {
                    normal_axis = i;
                }
            }
            if (normal_axis == -1)
            {
                throw _UTIL_MSG_EXCEPT(
                  "The boundary polygon did not align closely enough with a coordinate plane.");
            }

            // Find if we are closer to the negative or positive extent
            auto p             = element.node(0).point();
            double dist_to_neg = fabs(p[normal_axis] - extent.first[normal_axis]);
            double dist_to_pos = fabs(p[normal_axis] - extent.second[normal_axis]);
            double new_value;
            if (dist_to_neg < dist_to_pos)
            {
                new_value = extent.first[normal_axis];
            }
            else
            {
                new_value = extent.second[normal_axis];
            }
            // Move the nodes to lie on the plane
            for (auto& node_index : element.connectivity())
            {
                auto& node              = tow_surface.mesh.nodes()[node_index];
                node.point(normal_axis) = new_value;
            }
        }
    }
}



struct SegmentStats
{
    double mean;
    double std_dev;
};

SegmentStats calculate_segment_length_stats(VTMSTowMesh& tow_surface)
{
    std::vector<double> segment_lengths;
    for (auto& element : tow_surface.mesh.elements())
    {
        auto& conn = element.connectivity();
        for (auto segment : polygon_segments(element))
        {
            auto p1 = tow_surface.mesh.nodes()[conn[segment.first]].point();
            auto p2 = tow_surface.mesh.nodes()[conn[segment.second]].point();

            double segment_length = simple_math::distance(p1, p2);
            segment_lengths.push_back(segment_length);
        }
    }

    SegmentStats stats;
    stats.mean    = simple_math::calculate_mean(segment_lengths);
    stats.std_dev = simple_math::calculate_std_dev(segment_lengths);
    return stats;
};

SegmentStats calculate_segment_length_stats(const ElementConstAccessor& polygon)
{
    assert(polygon.type() == ElementType::POLYGON);

    std::vector<double> segment_lengths;
    auto& conn = polygon.connectivity();
    for (auto segment : polygon_segments(polygon))
    {
        auto p1 = polygon.node(conn[segment.first]).point();
        auto p2 = polygon.node(conn[segment.second]).point();

        double segment_length = simple_math::distance(p1, p2);
        segment_lengths.push_back(segment_length);
    }

    SegmentStats stats;
    stats.mean    = simple_math::calculate_mean(segment_lengths);
    stats.std_dev = simple_math::calculate_std_dev(segment_lengths);
    return stats;
};

/// Loops over the segments and removes any that are smaller than
/// factor*avg_segment_length and combines the nodes.
/// It then cleans up any polygons with two or less segments,
/// so there are no unrealistic polygons that are collapsed.
void remove_short_segments(VTMSTowMesh& tow_surface, const double factor)
{
    // Find the mean and standard deviation
    auto seg_stats    = calculate_segment_length_stats(tow_surface);
    double length_min = factor * seg_stats.mean;

    // Map, key: node to be replaced, value: node that replaces it
    std::map<LocalIndex, LocalIndex> node_replacements;

    {
        NodeToElementMap node_to_element;
        node_to_element.update(tow_surface.mesh);

        // Loop over the elements and remove edges that are too short
        for (auto& element : tow_surface.mesh.elements())
        {
            bool need_to_check_edges = true;
            while (need_to_check_edges)
            {
                // Set the flag to false and it will be set back to true if an edge is removed
                need_to_check_edges = false;

                // IF there is 2 or less points, then this element will be removed anyway.
                if (element.connectivity().size() < 3)
                {
                    break;
                }

                // Loop over the polygon segments and check if one is too short, which would need to
                // be removed
                auto& conn = element.connectivity();
                for (auto segment : polygon_segments(element))
                {
                    auto p1 = tow_surface.mesh.nodes()[conn[segment.first]].point();
                    auto p2 = tow_surface.mesh.nodes()[conn[segment.second]].point();

                    double segment_length = simple_math::distance(p1, p2);

                    if (segment_length < length_min)
                    {
                        // Which ever node lies on the fewest clipped planes will be removed and
                        // replace with the other
                        std::vector<Vector3D> boundary_normals_for_1;
                        for (auto element_index :
                             node_to_element.get_element_indices(conn[segment.first]))
                        {
                            if (tow_surface.surface_element_types[element_index] ==
                                VTMSTowMesh::SurfaceElementType::BOUNDARY)
                            {
                                auto normal = calculate_polygon_unit_normal(
                                  tow_surface.mesh.elements()[element_index]);
                                bool close_to_existing_normal = false;
                                for (auto& n : boundary_normals_for_1)
                                {
                                    if ((n - normal).isMuchSmallerThan(1.0, 1e-6))
                                    {
                                        close_to_existing_normal = true;
                                        break;
                                    }
                                }
                                if (!close_to_existing_normal)
                                {
                                    boundary_normals_for_1.push_back(normal);
                                }
                            }
                        }
                        std::vector<Vector3D> boundary_normals_for_2;
                        for (auto element_index :
                             node_to_element.get_element_indices(conn[segment.second]))
                        {
                            if (tow_surface.surface_element_types[element_index] ==
                                VTMSTowMesh::SurfaceElementType::BOUNDARY)
                            {
                                auto normal = calculate_polygon_unit_normal(
                                  tow_surface.mesh.elements()[element_index]);
                                bool close_to_existing_normal = false;
                                for (auto& n : boundary_normals_for_2)
                                {
                                    if ((n - normal).isMuchSmallerThan(1.0, 1e-6))
                                    {
                                        close_to_existing_normal = true;
                                        break;
                                    }
                                }
                                if (!close_to_existing_normal)
                                {
                                    boundary_normals_for_2.push_back(normal);
                                }
                            }
                        }

                        // If the size of the boundary normals are the same, then the algorithm has
                        // encountered a segment straddling two boundaries, such as segment 'a'
                        // connecting 1-3.  Skip this segment and let 1->2 and 3->2 later to
                        // collapse the element.
                        //
                        // (boundary)
                        // 1    2
                        // o----o
                        //  \   |
                        //   \  | (boundary)
                        //   a\ |
                        //     \|
                        //      o 3
                        if (boundary_normals_for_2.size() == boundary_normals_for_1.size())
                        {
                            continue;
                        }

                        LocalIndex node_to_remove;
                        LocalIndex node_to_keep;
                        if (boundary_normals_for_2.size() > boundary_normals_for_1.size())
                        {

                            node_to_keep   = conn[segment.second];
                            node_to_remove = conn[segment.first];
                        }
                        else
                        {
                            node_to_keep   = conn[segment.first];
                            node_to_remove = conn[segment.second];
                        }
                        node_replacements[node_to_remove] = node_to_keep;

                        // Remove the respective node from the connectivity
                        std::vector<LocalIndex> new_conn(conn.begin(), conn.end());
                        new_conn.erase(std::find(new_conn.begin(), new_conn.end(), node_to_remove));
                        tow_surface.mesh.elements().replace_connectivity(element.local_index(),
                                                                         utility::slice(new_conn));

                        // Loop over the elements that contain the node to remove
                        for (auto& parent_index :
                             node_to_element.get_element_indices(node_to_remove))
                        {
                            if (parent_index != element.local_index())
                            {
                                auto& parent = tow_surface.mesh.elements()[parent_index];

                                // If the element contains both nodes being combined, then delete
                                // the node being removed from the connectivity
                                auto& parent_conn = parent.connectivity();
                                if (parent_conn.contains(node_to_remove) &&
                                    parent_conn.contains(node_to_keep))
                                {
                                    std::vector<LocalIndex> new_parent_conn(parent_conn.begin(),
                                                                            parent_conn.end());
                                    new_parent_conn.erase(std::find(new_parent_conn.begin(),
                                                                    new_parent_conn.end(),
                                                                    node_to_remove));
                                    tow_surface.mesh.elements().replace_connectivity(
                                      parent.local_index(),
                                      utility::slice(new_parent_conn));
                                }
                                // If the element just contains node2, then just swap the nodes
                                else
                                {
                                    std::vector<LocalIndex> new_parent_conn(parent_conn.begin(),
                                                                            parent_conn.end());
                                    *std::find(new_parent_conn.begin(),
                                               new_parent_conn.end(),
                                               node_to_remove) = node_to_keep;
                                    tow_surface.mesh.elements().replace_connectivity(
                                      parent.local_index(),
                                      utility::slice(new_parent_conn));
                                }
                            }
                        }

                        // Debug
                        for (auto& element_2 : tow_surface.mesh.elements())
                        {
                            for (auto& i : element_2.connectivity())
                            {
                                assert(i != node_to_remove);
                            }
                        }

                        // Update parents so that node_to_element remains valid
                        node_to_element.adjust_for_node_combination(node_to_keep, node_to_remove);

                        // Break out of the loop over segments to get back into the while loop
                        // This will recreate the segments, which will be checked again.
                        need_to_check_edges = true;

                        break;

                    } // End if length < min length
                }     // End loop over edges
            }         // End while loop
        }             // End loop over elements
    }                 // End scope that removes segements

    // Remove the nodes that were merged into others
    LocalIndex nodes_removed_count = 0;
    std::vector<LocalIndex> old_to_new_node_index_map(tow_surface.mesh.nodes().size());
    std::iota(old_to_new_node_index_map.begin(), old_to_new_node_index_map.end(), 0);
    auto curr_map_index = 0;
    for (auto& replace_pair : node_replacements)
    {
        while (curr_map_index != replace_pair.first)
        {
            old_to_new_node_index_map[curr_map_index++] -= nodes_removed_count;
        }
        tow_surface.mesh.nodes().remove_stable(replace_pair.first - nodes_removed_count++);
        old_to_new_node_index_map[curr_map_index++] = replace_pair.second;
    }
    for (auto i = curr_map_index; i < old_to_new_node_index_map.size(); i++)
    {
        old_to_new_node_index_map[i] -= nodes_removed_count;
    }
    tow_surface.mesh.elements().permute_connectivity(old_to_new_node_index_map);

    // Remove the elements that only have 2 nodes (too many edges were collapsed)
    LocalIndex elements_removed_count = 0;
    for (auto element_index : utility::range(tow_surface.mesh.elements().size()))
    {
        auto& element = tow_surface.mesh.elements()[element_index - elements_removed_count];
        if (element.connectivity().size() <= 2)
        {
            // Account for the fact the local numbers are changing as elements are removed
            tow_surface.mesh.elements().remove_stable(element_index - elements_removed_count);
            tow_surface.surface_element_types.erase(tow_surface.surface_element_types.begin() +
                                                    (element_index - elements_removed_count));
            elements_removed_count++;
        }
    }

    // Check connectivities are at least within a plausible range (catching dumb mistakes)
    for (auto& element : tow_surface.mesh.elements())
    {
        for (auto& i : element.connectivity())
        {
            _UTIL_ASSERT(i >= 0 && i < tow_surface.mesh.nodes().size(),
                         "Invalid connectivity found for element " +
                           std::to_string(element.local_index()) + ", debug algorithm.");
        }
    }
};

/// Refines the segments that define the polygon, for segments larger
/// than the mean segment length plus numStdDevTol times the standard
/// deviation of the segment lengths.
void refine_polygon_segments(VTMSTowMesh& tow_surface, double length_max)
{
    // Loop over the polygons and refine the segments
    std::map<std::size_t, std::vector<LocalIndex>> refined_segments;
    std::map<std::size_t, LocalIndex> refinement_start_node;
    for (auto& element : tow_surface.mesh.elements())
    {
        // Loop over the edges and check if one is too long, which would need to be refined
        auto segments = polygon_segments(element);
        auto n_segments = segments.size();
        auto segment_index = 0;
        while (segment_index < n_segments)
        {
            auto& conn    = element.connectivity();
            auto& segment = segments[segment_index];
            auto p1       = tow_surface.mesh.nodes()[conn[segment.first]].point();
            auto p2       = tow_surface.mesh.nodes()[conn[segment.second]].point();

            double segment_length = simple_math::distance(p1, p2);

            if (segment_length > length_max)
            {
                // This hash will be the same regardless of the order
                auto hash = utility::hash_sorted(conn[segment.first], conn[segment.second]);

                // Check to see if this segment has already been refined
                auto find_iter = refined_segments.find(hash);
                if (find_iter != refined_segments.end())
                {
                    // Just add the nodes to the connectivity if it has been refined already
                    // Check to ensure the order is the correct direction
                    if (refinement_start_node[hash] == conn[segment.first])
                    {
                        std::vector<LocalIndex> new_conn(conn.begin(), conn.end());
                        new_conn.insert(
                          std::find(new_conn.begin(), new_conn.end(), conn[segment.first]) + 1,
                          find_iter->second.begin(),
                          find_iter->second.end());

                        tow_surface.mesh.elements().replace_connectivity(element.local_index(),
                                                                         utility::slice(new_conn));

                        // Update the segments since new points were inserted
                        segments = polygon_segments(element);
                        segment_index += find_iter->second.size() + 1;
                        n_segments = segments.size();
                    }
                    // Otherwise we need to reverse the order
                    else
                    {
                        auto reversed_nodes = find_iter->second;
                        std::reverse(reversed_nodes.begin(), reversed_nodes.end());

                        std::vector<LocalIndex> new_conn(conn.begin(), conn.end());
                        new_conn.insert(
                          std::find(new_conn.begin(), new_conn.end(), conn[segment.first]) + 1,
                          reversed_nodes.begin(),
                          reversed_nodes.end());

                        tow_surface.mesh.elements().replace_connectivity(element.local_index(),
                                                                         utility::slice(new_conn));

                        // Update the segments since new points were inserted
                        segments = polygon_segments(element);
                        segment_index += reversed_nodes.size() + 1;
                        n_segments = segments.size();
                    }
                }
                // If it hasn't been refined yet, then add nodes
                else
                {
                    int num_sub_segments =
                      static_cast<int>((segment_length + 0.5 * length_max) / length_max);
                    std::vector<Point> linspace_points =
                      simple_math::linspace(p1, p2, num_sub_segments);
                    std::vector<Point> new_points(linspace_points.begin() + 1,
                                                  linspace_points.end() - 1);
                    // Add the new nodes
                    std::vector<LocalIndex> new_nodes;
                    for (auto& new_point : new_points)
                    {
                        tow_surface.mesh.nodes().add(new_point);
                        new_nodes.push_back(tow_surface.mesh.nodes().size() - 1);
                    }

                    // Save the info for future segments
                    refined_segments[hash]      = new_nodes;
                    refinement_start_node[hash] = conn[segment.first];

                    // Add the new nodes to the element connectivity
                    std::vector<LocalIndex> new_conn(conn.begin(), conn.end());
                    new_conn.insert(
                      std::find(new_conn.begin(), new_conn.end(), conn[segment.first]) + 1,
                      new_nodes.begin(),
                      new_nodes.end());
                    tow_surface.mesh.elements().replace_connectivity(element.local_index(),
                                                                     utility::slice(new_conn));

                    // Update the segments since new points were inserted
                    segments = polygon_segments(element);
                    segment_index += num_sub_segments;
                    n_segments = segments.size();
                }
            } // End if length > max length
            else
            {
                segment_index++;
            }
        }     // End loop over edges
    }         // End loop over elements
} // End refine_polygon_segments

/// @brief TODO
/// @param output
/// @param points - list of points to include in the triangulation (assumes 2D, so 3rd coordinate is
///        ignored)
/// @param point_markers
/// @param boundary_segments
/// @param boundary_segment_markers
/// @param hole_points - list of points that indicating regions where a hole should lie (assumes 2D,
///        so 3rd coordinate is ignored).  Note: the point just has to be in any closed region
///        defined by the boundary segments.
/// @param attribute_points
/// @param attributes
/// @param allow_new_nodes_on_boundary
/// @param max_triangle_area
/// @param min_angle
triangulateio triangle_kernel(
  const std::vector<Point>& points,
  const std::vector<int>& point_markers,
  const std::vector<std::pair<LocalIndex, LocalIndex>>& boundary_segments,
  const std::vector<int>& boundary_segment_markers,
  const std::vector<std::vector<double>>& hole_points,
  const std::vector<std::vector<double>>& attribute_points,
  const std::vector<double>& attributes,
  const bool allow_new_nodes_on_boundary = false,
  const double max_triangle_area         = -1.0,
  double min_angle                       = 28.6)
{
    // Note - This method uses Triangle, a program by Jonathan Richard Shewchuk at Berkely
    // See his website at http://www.cs.cmu.edu/~quake/triangle.html for more information
    // This method is fashioned after tricall.c
    struct triangulateio in;

    in.numberofpoints          = static_cast<int>(points.size());
    in.numberofpointattributes = 0;
    in.pointlist               = new double[points.size() * 2];
    if (point_markers.size() > 0)
    {
        in.pointmarkerlist = new int[in.numberofpoints];
    }
    else
    {
        in.pointmarkerlist = NULL;
    }
    for (unsigned int i = 0; i < points.size(); ++i)
    {
        in.pointlist[i * 2]     = points[i][0];
        in.pointlist[i * 2 + 1] = points[i][1];
        if (point_markers.size() > 0)
        {
            in.pointmarkerlist[i] = point_markers[i];
        }
    }
    in.pointattributelist = NULL;
    in.numberofsegments   = static_cast<int>(boundary_segments.size());
    in.segmentlist        = new int[in.numberofsegments * 2];
    if (boundary_segment_markers.size() > 0)
    {
        in.segmentmarkerlist = new int[in.numberofsegments];
    }
    else
    {
        in.segmentmarkerlist = NULL;
    }
    for (int i = 0; i < in.numberofsegments; ++i)
    {
        in.segmentlist[i * 2]     = boundary_segments[i].first;
        in.segmentlist[i * 2 + 1] = boundary_segments[i].second;
        if (boundary_segment_markers.size() > 0)
        {
            in.segmentmarkerlist[i] = boundary_segment_markers[i];
        }
    }
    in.numberofholes = static_cast<int>(hole_points.size());
    in.holelist      = new double[hole_points.size() * 2];
    for (unsigned int i = 0; i < hole_points.size(); ++i)
    {
        in.holelist[i * 2]     = hole_points[i][0];
        in.holelist[i * 2 + 1] = hole_points[i][1];
    }
    in.numberofregions = static_cast<int>(attribute_points.size());
    if (attribute_points.size() != attributes.size())
    {
        delete[] in.pointlist;
        delete[] in.segmentlist;
        delete[] in.holelist;
        throw _UTIL_MSG_EXCEPT("AttributeCoords and Attributes are not same size!");
    }
    in.regionlist = new double[attribute_points.size() * 4];
    for (unsigned int i = 0; i < attribute_points.size(); ++i)
    {
        in.regionlist[i * 4]     = attribute_points[i][0];
        in.regionlist[i * 4 + 1] = attribute_points[i][1];
        in.regionlist[i * 4 + 2] = attributes[i];
        in.regionlist[i * 4 + 3] = -1.0; // Area constraint, negative means none applied.
    }

    triangulateio output;
    output.pointlist             = NULL;
    output.pointattributelist    = NULL;
    output.pointmarkerlist       = NULL;
    output.trianglelist          = NULL;
    output.triangleattributelist = NULL;
    output.neighborlist          = NULL;
    output.segmentlist           = NULL;
    output.segmentmarkerlist     = NULL;
    output.edgelist              = NULL;
    output.edgemarkerlist        = NULL;

    std::ostringstream tri_arg_stream;
    tri_arg_stream << "AQp";
    // Note: a double Y indicates that Triangle will not add any new points.
    // One Y means prohibits the any new nodes on the boundary
    if (!allow_new_nodes_on_boundary)
    {
        tri_arg_stream << "Y";
    }
    tri_arg_stream << "ezDq" << min_angle;
    if (max_triangle_area > 0)
    {
        tri_arg_stream.precision(15);
        tri_arg_stream << "a" << std::fixed << max_triangle_area;
    }
    char args[50];
    strcpy_s(args, tri_arg_stream.str().c_str());

    triangulate(args, &in, &output, (struct triangulateio*)NULL);

    // frees memory used by triangle
    delete[] in.pointlist;
    delete[] in.segmentlist;
    if (in.segmentmarkerlist != NULL)
    {
        delete[] in.segmentmarkerlist;
    }
    delete[] in.holelist;
    delete[] in.regionlist;
    free(in.pointattributelist);

    return output;
}

void triangulate_polygon(VTMSTowMesh& tow_surface,
                         const ElementConstAccessor& polygon_element,
                         const double max_area)
{
    auto n_nodes = polygon_element.connectivity().size();

    // Define a coordinate system
    auto mean_point = calc_mean_coordinate(polygon_element);
    auto p1         = polygon_element.node(0).point();
    auto p2         = polygon_element.node(1).point();
    auto axis_1     = simple_math::point_as_vector(p1 - mean_point).normalized();
    auto axis_3     = calculate_polygon_unit_normal(polygon_element);
    auto axis_2     = axis_3.cross(axis_1);

    // Transform coordinates to this normalized coordinate system
    simple_math::Matrix3D transformation;
    transformation << axis_1[0], axis_1[1], axis_1[2], axis_2[0], axis_2[1], axis_2[2], axis_3[0],
      axis_3[1], axis_3[2];
    std::vector<Point> transformed_points;
    transformed_points.reserve(polygon_element.connectivity().size());
    // Save the nodes' indices in the connectivity as markers so we can recover the nodes after
    // triangulation
    std::vector<int> markers;
    markers.reserve(polygon_element.connectivity().size());
    double normalized_z_offset = 0.0;
    LocalIndex count           = 0;
    for (auto& node : polygon_element.nodes())
    {
        auto transformed_p  = transformation * node.coord();
        normalized_z_offset = transformed_p(2);
        transformed_points.push_back({ transformed_p(0), transformed_p(1), 0.0 });
        // n->set_local_number(count);
        // uniqueNodesVec.push_back(n);
        markers.push_back(count);
        count++;
    }

    // Create the boundary segments for the triangulation
    std::vector<std::pair<LocalIndex, LocalIndex>> boundary_segments;
    for (auto segment : polygon_segments(polygon_element))
    {
        boundary_segments.push_back({ segment.first, segment.second });
    }

    // Call triangle
    struct triangulateio output = triangle_kernel(transformed_points,
                                                  markers,
                                                  boundary_segments,
                                                  {},
                                                  {},
                                                  {},
                                                  {},
                                                  false,
                                                  max_area,
                                                  10.0);

    {
        // So, I will exploit the fact that triangle labels the given points first
        // Add all the points, starting after the points I seeded
        transformation.transposeInPlace(); // Now this transforms normalized to global
        std::vector<LocalIndex> nodes_in_triangulation(polygon_element.connectivity().begin(),
                                                       polygon_element.connectivity().end());
        for (int i = static_cast<int>(polygon_element.connectivity().size());
             i < output.numberofpoints;
             i++)
        {
            // Note: we need to transform the coords back to global system
            simple_math::Vector3D coord_2D = { output.pointlist[i * 2],
                                               output.pointlist[i * 2 + 1],
                                               normalized_z_offset };
            auto global_p                  = simple_math::to_point(transformation * coord_2D);
            tow_surface.mesh.nodes().add(global_p);
            nodes_in_triangulation.push_back(tow_surface.mesh.nodes().size() - 1);
        }
        // Now form the connectivities for the new triangles
        std::vector<LocalIndex> conn(3, 0);
        for (int i = 0; i < output.numberoftriangles; i++)
        {
            for (auto j : utility::range(3))
            {
                conn[j] = nodes_in_triangulation[output.trianglelist[i * 3 + j]];
            }
            tow_surface.mesh.elements().add(ElementType::T2D3N, conn);
        }
    }

    // frees memory used by triangle
    free(output.pointlist);
    free(output.pointattributelist);
    free(output.pointmarkerlist);
    free(output.trianglelist);
    free(output.triangleattributelist);
    free(output.neighborlist);
    free(output.segmentlist);
    free(output.segmentmarkerlist);
    free(output.edgelist);
    free(output.edgemarkerlist);
}

void submesh_polygons(VTMSTowMesh& tow_surface, double max_tri_area)
{
    // Loop over elements in surface mesh and submesh polygons with more than 4 nodes
    LocalIndex n_elements_removed = 0;
    for (auto polygon_element_index : utility::range(tow_surface.mesh.elements().size()))
    {
        auto& element = tow_surface.mesh.elements()[polygon_element_index - n_elements_removed];
        if (element.connectivity().size() > 3)
        {
            std::vector<std::vector<LocalIndex>> boundary_segments;
            auto& conn = element.connectivity();
            for (auto segment : polygon_segments(element))
            {
                boundary_segments.push_back({ conn[segment.first], conn[segment.second] });
            }
            triangulate_polygon(tow_surface, element, max_tri_area);

            // Remove the original polygon element
            tow_surface.mesh.elements().remove_stable(polygon_element_index - n_elements_removed++);
        }
        else if (element.connectivity().size() == 3)
        {
            // Convert polygons to trianlges
            element.type() = ElementType::T2D3N;
        }
        else
        {
            throw _UTIL_MSG_EXCEPT("Bad element in surface.");
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Definitions from header

VTMSTowMesh read_vtms_ctw(const std::string& file_path)
{
    // Open the ctw file
    std::ifstream fin(file_path.c_str());
    if (!fin)
    {
        throw _UTIL_MSG_EXCEPT("Cannot find file " + file_path);
    }

    VTMSTowMesh tow_surface;
    tow_surface.type = VTMSTowMesh::MeshType::POLY;

    // Read the number of nodes and material number
    std::string token;
    fin >> token >> token >> token;
    fin >> tow_surface.material_number;
    fin >> token;
    int n_nodes;
    fin >> n_nodes;

    // Read the coords from polygon file, but these will be added to as polygons are submeshed
    int number;
    double x, y, z;
    for (int i = 0; i < n_nodes; i++)
    {
        fin >> number >> x >> y >> z;
        tow_surface.mesh.nodes().add({ x, y, z });
    }

    // Read the polygons
    int n_polys;
    fin >> token >> n_polys;
    int n_vertices;
    std::string line;
    tow_surface.surface_element_types.reserve(n_polys);
    for (int i = 0; i < n_polys; i++)
    {
        // VTMS is not consistent with the text format and adds a node
        // set sometimes.  This is a workaround.
        std::getline(fin, line);
        while (line.find_first_of("0123456789") == std::string::npos)
        {
            std::getline(fin, line);
        }

        auto tokens = utility::tokenize(line);
        number      = atoi(tokens[0].c_str());
        n_vertices  = atoi(tokens[1].c_str());

        LocalIndex node_index;
        std::vector<LocalIndex> connectivity(n_vertices, 0);
        for (int j = 0; j < n_vertices; j++)
        {
            fin >> node_index;
            connectivity[j] = node_index;
        }

        // Note: this is a hack to determine if the polygonal is on the boundary or not.  Ideally,
        // VTMS should output this information, but for now, we will assume if it has more than 10
        // points, then it was clipped.
        if (connectivity.size() > 10)
        {
            tow_surface.surface_element_types.push_back(VTMSTowMesh::SurfaceElementType::BOUNDARY);
        }
        else
        {
            tow_surface.surface_element_types.push_back(VTMSTowMesh::SurfaceElementType::INTERIOR);
        }

        tow_surface.mesh.elements().add(ElementType::POLYGON, connectivity);
    }

    // It will start pointing outward, but we will switch this during algorithms
    tow_surface.surface_normal_direction = VTMSTowMesh::NormalDirection::OUTWARD;

    // Skip to the medials data
    line = utility::skip_to_line_that_contains("*Medials:", fin, line);

    // Now we are at the medials data
    unsigned int n_medials;
    fin >> n_medials;

    // We cannot preallocate since VTMS sometimes gives indefinite vlaues
    unsigned int num;
    std::string val;
    for (unsigned int i = 0; i < n_medials; i++)
    {
        fin >> num;
        Point p;
        bool found_indefinite = false;
        for (unsigned int j = 0; j < 3; j++)
        {
            fin >> val;
            if (utility::istr_is_equal(val, "-1.#IND00000") == false)
            {
                p[j] = atof(val.c_str());
            }
            else
            {
                found_indefinite = true;
            }
        }

        if (!found_indefinite)
        {
            tow_surface.medials.push_back(p);
        }
    }

    snap_boundary_points(tow_surface);

    return tow_surface;
}

double triangulate_tow_surface(VTMSTowMesh& tow_surface)
{
    _UTIL_ASSERT(tow_surface.type == VTMSTowMesh::MeshType::POLY,
                 "The tow surface must be polygonal.");
    _UTIL_ASSERT(tow_surface.surface_element_types.size() == tow_surface.mesh.elements().size(),
                 "The surface element types must match the number of elements in the mesh.");

    // Remove really short segments that would result in bad elements
    // Will remove 25% of the average segment length or smaller
    remove_short_segments(tow_surface, 0.25);

    // write_vtk(tow_surface.mesh, "surface_post_remove_short_segments.vtk");

    // Find the mean and standard deviation
    auto n_std_devs   = 4;
    auto seg_stats    = calculate_segment_length_stats(tow_surface);
    double length_max = seg_stats.mean + n_std_devs * seg_stats.std_dev;

    // Refine to ensure segments are within two std devs of the mean
    refine_polygon_segments(tow_surface, length_max);

    // write_vtk(tow_surface.mesh, "surface_post_refine_polygon_segments.vtk");

    seg_stats = calculate_segment_length_stats(tow_surface);
    submesh_polygons(tow_surface, seg_stats.mean * seg_stats.mean);

    return seg_stats.mean;
}

void write_to_ply(const Mesh& triangular_surface, const std::string& file_path)
{
    std::ofstream fout(file_path.c_str());
    if (!fout)
    {
        throw _UTIL_MSG_EXCEPT("Could not create file " + file_path);
    }

    // Write header section, refer to http://paulbourke.net/dataformats/ply/
    fout << "ply" << std::endl;
    fout << "format ascii 1.0" << std::endl;
    fout << "comment Mesh exported from mesh module of SHEFF framework" << std::endl;
    fout << "element vertex " << triangular_surface.nodes().size() << std::endl;
    fout << "property float x" << std::endl;
    fout << "property float y" << std::endl;
    fout << "property float z" << std::endl;
    fout << "element face " << triangular_surface.elements().size() << std::endl;
    fout << "property list uchar int vertex_index" << std::endl;
    fout << "end_header" << std::endl;

    // Write vertex list
    for (auto& node : triangular_surface.nodes())
    {
        fout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << node.x() << " "
             << node.y() << " " << node.z() << std::endl;
    }

    // Write face list
    for (auto& element : triangular_surface.elements())
    {
        _UTIL_ASSERT(element.type() == ElementType::T2D3N,
                     "The tow surface must only contain triangular elements.");

        fout << "3 ";
        for (auto& i : element.connectivity())
        {
            fout << i << " ";
        }
        fout << std::endl;
    }
    fout.close();
}

Mesh create_triangular_matrix_box(const Point& n,
                                  const Point& p,
                                  const double& target_segment_length)
{
    //     ^z
    //
    //     4---7
    //    /|  /|
    //   5-0-6 3 --> y
    //   |/  |/
    //   1---2
    //
    //  / x
    //
    // where 0 is most negative point, and 6 is most positive point
    Mesh box;

    box.nodes().add(n);
    box.nodes().add({ p[0], n[1], n[2] });
    box.nodes().add({ p[0], p[1], n[2] });
    box.nodes().add({ n[0], p[1], n[2] });
    box.nodes().add({ n[0], n[1], p[2] });
    box.nodes().add({ p[0], n[1], p[2] });
    box.nodes().add(p);
    box.nodes().add({ n[0], p[1], p[2] });

    if (target_segment_length > 0.0)
    {
        box.elements().add(ElementType::POLYGON, { 0, 3, 2, 1 }); // Bottom
        box.elements().add(ElementType::POLYGON, { 1, 2, 6, 5 }); // Front
        box.elements().add(ElementType::POLYGON, { 0, 4, 7, 3 }); // Back
        box.elements().add(ElementType::POLYGON, { 1, 5, 4, 0 }); // Left side
        box.elements().add(ElementType::POLYGON, { 3, 7, 6, 2 }); // Right side
        box.elements().add(ElementType::POLYGON, { 6, 7, 4, 5 }); // Top

        VTMSTowMesh tmp_vtms_obj;
        tmp_vtms_obj.mesh = std::move(box);
        for (auto& element : tmp_vtms_obj.mesh.elements())
        {
            tmp_vtms_obj.surface_element_types.push_back(VTMSTowMesh::SurfaceElementType::BOUNDARY);
        }
        tmp_vtms_obj.type = VTMSTowMesh::MeshType::POLY;

        // Refine to ensure segments are within two std devs of the mean
        refine_polygon_segments(tmp_vtms_obj, target_segment_length);

        //write_vtk(tmp_vtms_obj.mesh, "surface_post_refine_polygon_segments.vtk");

        submesh_polygons(tmp_vtms_obj, 0.5 * target_segment_length * target_segment_length);

        box = std::move(tmp_vtms_obj.mesh);
        return box;
    }
    else
    {
        box.elements().add(ElementType::T2D3N, { 2, 1, 0 }); // Bottom
        box.elements().add(ElementType::T2D3N, { 2, 0, 3 });
        box.elements().add(ElementType::T2D3N, { 6, 5, 1 }); // Front
        box.elements().add(ElementType::T2D3N, { 6, 1, 2 });
        box.elements().add(ElementType::T2D3N, { 0, 4, 7 }); // Back
        box.elements().add(ElementType::T2D3N, { 0, 7, 3 });
        box.elements().add(ElementType::T2D3N, { 1, 5, 0 }); // Left side
        box.elements().add(ElementType::T2D3N, { 0, 5, 4 });
        box.elements().add(ElementType::T2D3N, { 3, 7, 6 }); // Right side
        box.elements().add(ElementType::T2D3N, { 3, 6, 2 });
        box.elements().add(ElementType::T2D3N, { 6, 7, 4 }); // Top
        box.elements().add(ElementType::T2D3N, { 6, 4, 5 });

        return box;
    }
}
}