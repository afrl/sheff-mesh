#include "mesh/algorithms/Polygons.h"

#include "mesh/algorithms/IO.h"
#include "mesh/algorithms/Manipulation.h"
#include "mesh/kd_tree/KDTree.h"
#include "mesh/relations/NodeToElementMap.h"
#include "second_party_mesh/triangle/triangle.h"
#include "simple_math/Utility.h"
#include "utility/Hash.h"

#include <algorithm>
#include <string.h>

using namespace simple_math::point_operators;

namespace mesh
{
//////////////////////////////////////////////////////////////////////////////////////////////////
// Helpers

namespace detail
{

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
triangulateio triangle_polygon_kernel(
  const std::vector<Point>& points,
  const std::vector<int>& point_markers,
  const std::vector<std::pair<LocalIndex, LocalIndex>>& boundary_segments,
  const std::vector<int>& boundary_segment_markers,
  const std::vector<std::vector<double>>& hole_points,
  const std::vector<std::vector<double>>& attribute_points,
  const std::vector<double>& attributes,
  const bool allow_new_nodes_on_boundary,
  const double max_triangle_area,
  double min_angle)
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
    if (!allow_new_nodes_on_boundary && max_triangle_area < 0.0)
    {
        tri_arg_stream << "YY";
    }
    else if (!allow_new_nodes_on_boundary)
    {
        tri_arg_stream << "Y";
    }
    tri_arg_stream << "ezD";
    if (min_angle > 0.0 && max_triangle_area > 0.0)
    {
        tri_arg_stream << "q" << min_angle;
    }
    if (max_triangle_area > 0.0)
    {
        tri_arg_stream.precision(15);
        tri_arg_stream << "a" << std::fixed << max_triangle_area;
    }
    char args[50];
    strcpy(args, tri_arg_stream.str().c_str());

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
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Definitions from header

std::vector<std::pair<ConnIndex, ConnIndex>> polygon_segments(
  const utility::Slice<const LocalIndex, ConnIndex>& connectivity)
{
    std::vector<std::pair<ConnIndex, ConnIndex>> segments;
    for (auto i = 0; i < connectivity.size() - 1; i++)
    {
        segments.push_back({ i, i + 1 });
    }
    segments.push_back({ static_cast<ConnIndex>(connectivity.size() - 1), 0 });
    return segments;
}

std::vector<std::pair<ConnIndex, ConnIndex>> polygon_segments(
  const utility::Slice<LocalIndex, ConnIndex>& connectivity)
{
    return polygon_segments(
      utility::Slice<const LocalIndex, ConnIndex>(connectivity.data(), connectivity.size()));
}


// Implements Newell's method to calculate a polygons normal vector
simple_math::Vector3D calc_polygon_unit_normal(const ElementConstAccessor& polygon_element)
{
    simple_math::Vector3D normal = simple_math::Vector3D::Zero();
    for (auto segment : polygon_segments(polygon_element.connectivity()))
    {
        auto p0 = polygon_element.node(segment.first).coord();
        auto p1 = polygon_element.node(segment.second).coord();
        normal += p0.cross(p1);
    }
    return normal.normalized();
}

simple_math::Vector3D calc_polygon_unit_normal(const std::vector<simple_math::Point>& points)
{
    simple_math::Vector3D normal = simple_math::Vector3D::Zero();
    for (auto i : utility::range(points.size() - 1))
    {
        auto p0 = simple_math::point_as_vector(points[i]);
        auto p1 = simple_math::point_as_vector(points[i + 1]);
        normal += p0.cross(p1);
    }
    { // Closing the loop
        auto p0 = simple_math::point_as_vector(points.back());
        auto p1 = simple_math::point_as_vector(points.front());
        normal += p0.cross(p1);
    }
    return normal.normalized();
}

SegmentStats calc_segment_length_stats(const Mesh& m)
{
    std::vector<double> segment_lengths;
    for (auto element : m.elements())
    {
        auto conn = element.connectivity();
        for (auto segment : polygon_segments(element.connectivity()))
        {
            auto p1 = m.nodes()[conn[segment.first]].point();
            auto p2 = m.nodes()[conn[segment.second]].point();

            double segment_length = simple_math::distance(p1, p2);
            segment_lengths.push_back(segment_length);
        }
    }

    SegmentStats stats;
    stats.mean    = simple_math::calculate_mean(segment_lengths);
    stats.std_dev = simple_math::calculate_std_dev(segment_lengths);
    return stats;
}

SegmentStats calc_segment_length_stats(const ElementConstAccessor& polygon_element)
{
    assert(polygon_element.type() == ElementType::POLYGON);

    std::vector<double> segment_lengths;
    auto conn = polygon_element.connectivity();
    for (auto segment : polygon_segments(conn))
    {
        auto p1 = polygon_element.node(segment.first).point();
        auto p2 = polygon_element.node(segment.second).point();

        double segment_length = simple_math::distance(p1, p2);
        segment_lengths.push_back(segment_length);
    }

    SegmentStats stats;
    stats.mean    = simple_math::calculate_mean(segment_lengths);
    stats.std_dev = simple_math::calculate_std_dev(segment_lengths);
    return stats;
}

CollapseSegmentsResult collapse_short_polygon_segments(
  Mesh& mesh,
  const double min_edge_length,
  const std::vector<std::pair<LocalIndex, LocalIndex>>& feature_edges,
  const std::vector<LocalIndex>& edge_id_to_feature_ids_map)
{
    // Map, key: node to be replaced, value: node that replaces it
    std::map<LocalIndex, LocalIndex> node_replacements;

    // A map from the node index to (pair first) the feature ids of all feature edges that contain
    // the node and (pair second) the other node in the feature edge.
    std::unordered_map<LocalIndex, std::vector<std::pair<LocalIndex, LocalIndex>>>
      node_to_feature_id_and_other_node_id_map;
    for (auto i : utility::range(feature_edges.size()))
    {
        auto& feature_id = edge_id_to_feature_ids_map[i];
        auto& n1_id      = feature_edges[i].first;
        auto& n2_id      = feature_edges[i].second;

        auto find_iter1 = node_to_feature_id_and_other_node_id_map.find(n1_id);
        if (find_iter1 == node_to_feature_id_and_other_node_id_map.end())
        {
            node_to_feature_id_and_other_node_id_map[n1_id] = { { feature_id, n2_id } };
        }
        else
        {
            find_iter1->second.push_back({ feature_id, n2_id });
        }

        auto find_iter2 = node_to_feature_id_and_other_node_id_map.find(n2_id);
        if (find_iter2 == node_to_feature_id_and_other_node_id_map.end())
        {
            node_to_feature_id_and_other_node_id_map[n2_id] = { { feature_id, n1_id } };
        }
        else
        {
            find_iter2->second.push_back({ feature_id, n1_id });
        }
    }

    NodeToElementMap node_to_element;
    node_to_element.update(mesh);

    // Loop over the elements and remove edges that are too short
    for (auto element : mesh.elements())
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
            auto conn = element.connectivity();
            for (auto segment : polygon_segments(element.connectivity()))
            {
                const auto& n1_id = conn[segment.first];
                const auto& n2_id = conn[segment.second];


                auto p1 = mesh.nodes()[n1_id].point();
                auto p2 = mesh.nodes()[n2_id].point();

                double segment_length = simple_math::distance(p1, p2);

                if (segment_length < min_edge_length)
                {

                    LocalIndex node_to_remove;
                    LocalIndex node_to_keep;
                    auto find_iter_n1 = node_to_feature_id_and_other_node_id_map.find(n1_id);
                    auto find_iter_n2 = node_to_feature_id_and_other_node_id_map.find(n2_id);
                    if (find_iter_n1 != node_to_feature_id_and_other_node_id_map.end() &&
                        find_iter_n2 == node_to_feature_id_and_other_node_id_map.end())
                    { // case where n1 is part of a feature but n2 is not
                        node_to_remove = n2_id;
                        node_to_keep   = n1_id;
                    }
                    else if (find_iter_n1 == node_to_feature_id_and_other_node_id_map.end() &&
                             find_iter_n2 != node_to_feature_id_and_other_node_id_map.end())
                    { // case where n2 is part of a feature but n1 is not
                        node_to_remove = n1_id;
                        node_to_keep   = n2_id;
                    }
                    else if (find_iter_n1 == node_to_feature_id_and_other_node_id_map.end() &&
                             find_iter_n2 == node_to_feature_id_and_other_node_id_map.end())
                    { // neither is part of a feature
                        // It is arbitrary which one to remove
                        node_to_remove = n2_id;
                        node_to_keep   = n1_id;
                    }
                    else
                    { // Both are part of features, so the collapse needs to be done carefully
                        std::vector<LocalIndex> n1_features;
                        for (auto& v : find_iter_n1->second)
                        {
                            n1_features.push_back(v.first);
                        }
                        std::sort(n1_features.begin(), n1_features.end());
                        n1_features.erase(std::unique(n1_features.begin(), n1_features.end()),
                                          n1_features.end());

                        std::vector<LocalIndex> n2_features;
                        for (auto& v : find_iter_n2->second)
                        {
                            n2_features.push_back(v.first);
                        }
                        std::sort(n2_features.begin(), n2_features.end());
                        n2_features.erase(std::unique(n2_features.begin(), n2_features.end()),
                                          n2_features.end());

                        std::vector<LocalIndex> intersection;
                        std::set_intersection(n1_features.begin(),
                                              n1_features.end(),
                                              n2_features.begin(),
                                              n2_features.end(),
                                              std::back_inserter(intersection));
                        if (n1_features == intersection && n2_features == intersection)
                        {
                            // It is arbitrary which one to remove
                            node_to_remove = n2_id;
                            node_to_keep   = n1_id;
                        }
                        else if (n1_features == intersection && n2_features != intersection)
                        {
                            node_to_remove = n1_id;
                            node_to_keep   = n2_id;
                        }
                        else if (n1_features != intersection && n2_features == intersection)
                        {
                            node_to_remove = n2_id;
                            node_to_keep   = n1_id;
                        }
                        else
                        {
                            // The algorithm hasencountered a segment straddling two boundaries,
                            // such as segment 'a' connecting 1-3.  Skip this segment and let 1->2
                            // and 3->2 later to collapse the element.
                            //
                            // (feature 0)
                            // 1    2
                            // o----o
                            //  \   |
                            //   \  | (feature 1)
                            //   a\ |
                            //     \|
                            //      o 3
                            continue;
                        }
                    }

                    node_replacements[node_to_remove] = node_to_keep;

                    // Remove the respective node from the connectivity
                    std::vector<LocalIndex> new_conn(conn.begin(), conn.end());
                    new_conn.erase(std::find(new_conn.begin(), new_conn.end(), node_to_remove));
                    mesh.elements().replace_connectivity(element.local_index(),
                                                         utility::slice(new_conn));

                    // Loop over the elements that contain the node to remove
                    for (auto parent_index : node_to_element.get_element_indices(node_to_remove))
                    {
                        if (parent_index != element.local_index())
                        {
                            auto parent = mesh.elements()[parent_index];

                            // If the element contains both nodes being combined, then delete
                            // the node being removed from the connectivity
                            auto parent_conn = parent.connectivity();
                            if (parent_conn.contains(node_to_remove) &&
                                parent_conn.contains(node_to_keep))
                            {
                                std::vector<LocalIndex> new_parent_conn(parent_conn.begin(),
                                                                        parent_conn.end());
                                new_parent_conn.erase(std::find(new_parent_conn.begin(),
                                                                new_parent_conn.end(),
                                                                node_to_remove));
                                mesh.elements().replace_connectivity(
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
                                mesh.elements().replace_connectivity(
                                  parent.local_index(),
                                  utility::slice(new_parent_conn));
                            }
                        }
                    }

                    // Debug
                    for (auto element_2 : mesh.elements())
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

    // Remove the nodes that were merged into others
    std::vector<std::pair<LocalIndex, LocalIndex>> combinations;
    for (auto& combination : node_replacements)
    {
        combinations.push_back({ combination.first, combination.second });
    }
    combine_nodes(mesh, combinations);

    // Remove the elements that only have 2 nodes (too many edges were collapsed)
    LocalIndex elements_removed_count = 0;
    std::vector<LocalIndex> elements_removed;
    for (auto element_index : utility::range(mesh.elements().size()))
    {
        auto element = mesh.elements()[element_index - elements_removed_count];
        if (element.connectivity().size() <= 2)
        {
            // Account for the fact the local numbers are changing as elements are removed
            mesh.elements().remove_stable(element_index - elements_removed_count);
            elements_removed.push_back(element_index);
            elements_removed_count++;
        }
    }

    // Check connectivities are at least within a plausible range (catching dumb mistakes)
    for (auto element : mesh.elements())
    {
        for (auto& i : element.connectivity())
        {
            _UTIL_ASSERT(i >= 0 && i < mesh.nodes().size(),
                         "Invalid connectivity found for element " +
                           std::to_string(element.local_index()) + ", debug algorithm.");
        }
    }

    return { combinations, elements_removed };
}

std::unordered_map<std::size_t, std::vector<LocalIndex>> refine_all_polygon_segments(
  Mesh& mesh,
  double length_max)
{
    // Loop over the polygons and refine the segments
    std::unordered_map<std::size_t, std::vector<LocalIndex>> edge_hash_to_new_conn_map;
    for (auto& element : mesh.elements())
    {
        // Loop over the edges and check if one is too long, which would need to be refined
        auto segments      = polygon_segments(element.connectivity());
        auto n_segments    = segments.size();
        auto segment_index = 0;
        while (segment_index < n_segments)
        {
            auto conn     = element.connectivity();
            auto& segment = segments[segment_index];
            auto p1       = mesh.nodes()[conn[segment.first]].point();
            auto p2       = mesh.nodes()[conn[segment.second]].point();

            double segment_length = simple_math::distance(p1, p2);

            if (segment_length > length_max)
            {
                // This hash will be the same regardless of the order
                auto hash = utility::hash_sorted(conn[segment.first], conn[segment.second]);

                // Check to see if this segment has already been refined
                auto find_iter = edge_hash_to_new_conn_map.find(hash);
                if (find_iter != edge_hash_to_new_conn_map.end())
                {
                    std::vector<LocalIndex> new_polygon_conn(conn.begin(), conn.end());
                    // Just add the nodes to the connectivity if it has been refined already
                    // Check to ensure the order is the correct direction
                    if (find_iter->second[0] == conn[segment.first])
                    {
                        new_polygon_conn.insert(std::find(new_polygon_conn.begin(),
                                                          new_polygon_conn.end(),
                                                          conn[segment.first]) +
                                                  1,
                                                find_iter->second.begin() + 1,
                                                find_iter->second.end() - 1);

                        mesh.elements().replace_connectivity(element.local_index(),
                                                             utility::slice(new_polygon_conn));
                    }
                    // Otherwise we need to reverse the order
                    else
                    {
                        auto reversed_nodes = find_iter->second;
                        std::reverse(reversed_nodes.begin(), reversed_nodes.end());

                        new_polygon_conn.insert(std::find(new_polygon_conn.begin(),
                                                          new_polygon_conn.end(),
                                                          conn[segment.first]) +
                                                  1,
                                                reversed_nodes.begin() + 1,
                                                reversed_nodes.end() - 1);

                        mesh.elements().replace_connectivity(element.local_index(),
                                                             utility::slice(new_polygon_conn));
                    }

                    // Update the segments since new points were inserted
                    segments = polygon_segments(element.connectivity());
                    segment_index += static_cast<int>(find_iter->second.size()) - 1;
                    n_segments = segments.size();
                }
                // If it hasn't been refined yet, then add nodes
                else
                {
                    int num_sub_segments =
                      static_cast<int>((segment_length + 0.5 * length_max) / length_max);
                    if (num_sub_segments < 2)
                    {
                        num_sub_segments = 2;
                    }
                    std::vector<Point> linspace_points =
                      simple_math::linspace(p1, p2, num_sub_segments);
                    std::vector<Point> new_points(linspace_points.begin() + 1,
                                                  linspace_points.end() - 1);
                    // Add the new nodes
                    std::vector<LocalIndex> new_nodes;
                    for (auto& new_point : new_points)
                    {
                        mesh.nodes().add(new_point);
                        new_nodes.push_back(mesh.nodes().size() - 1);
                    }


                    // Save the info for future segments
                    std::vector<LocalIndex> new_edge_conn = { conn[segment.first],
                                                              conn[segment.second] };
                    new_edge_conn.insert(new_edge_conn.begin() + 1,
                                         new_nodes.begin(),
                                         new_nodes.end());
                    edge_hash_to_new_conn_map[hash] = new_edge_conn;

                    // Add the new nodes to the element connectivity
                    std::vector<LocalIndex> new_polygon_conn(conn.begin(), conn.end());
                    new_polygon_conn.insert(std::find(new_polygon_conn.begin(),
                                                      new_polygon_conn.end(),
                                                      conn[segment.first]) +
                                              1,
                                            new_nodes.begin(),
                                            new_nodes.end());
                    mesh.elements().replace_connectivity(element.local_index(),
                                                         utility::slice(new_polygon_conn));

                    // Update the segments since new points were inserted
                    segments = polygon_segments(element.connectivity());
                    segment_index += num_sub_segments;
                    n_segments = segments.size();
                }
            } // End if length > max length
            else
            {
                segment_index++;
            }
        } // End loop over edges
    }     // End loop over elements
    return edge_hash_to_new_conn_map;
} // End refine_all_polygons

std::vector<LocalIndex> triangulate_polygon(Mesh& total_mesh,
                                            const ElementConstAccessor& polygon_element,
                                            const double max_area_a,
                                            const bool remove_polygon)
{
    assert(&polygon_element.node_container() == &total_mesh.nodes());

    double max_area      = max_area_a;
    auto orig_n_elements = total_mesh.elements().size();

    // Define a coordinate system
    auto mean_point = calc_mean_coordinate(polygon_element);
    auto p1         = polygon_element.node(0).point();
    auto p2         = polygon_element.node(1).point();
    auto axis_1     = simple_math::point_as_vector(p1 - mean_point).normalized();
    auto axis_3     = calc_polygon_unit_normal(polygon_element);
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
        markers.push_back(count);
        count++;
    }

    // Create the boundary segments for the triangulation
    std::vector<std::pair<LocalIndex, LocalIndex>> boundary_segments;
    for (auto segment : polygon_segments(polygon_element.connectivity()))
    {
        boundary_segments.push_back({ segment.first, segment.second });
    }

    if (max_area == std::numeric_limits<double>::max())
    {
        auto stats = calc_segment_length_stats(polygon_element);
        max_area   = stats.mean * stats.mean;
    }

    // Call triangle
    triangulateio output;
    try
    {
        output = detail::triangle_polygon_kernel(transformed_points,
                                                 markers,
                                                 boundary_segments,
                                                 {},
                                                 {},
                                                 {},
                                                 {},
                                                 false,
                                                 max_area,
                                                 10.0);
    }
    catch (...)
    {
        Mesh tmp;
        tmp.copy(total_mesh, polygon_element.global_index());

        for (auto& p : transformed_points)
        {
            tmp.nodes().add(p);
        }
        auto offset = polygon_element.connectivity().size();
        for (auto& s : boundary_segments)
        {
            tmp.elements().add(ElementType::L1D2N, { s.first + offset, s.second + offset });
        }
        write_vtk(tmp, "failed_triangle_input.vtk", false);
    }

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

            auto global_coord = transformation * coord_2D;
            if (!global_coord.array().isFinite().all())
            {
                Mesh tmp;
                tmp.copy(total_mesh, polygon_element.global_index());

                for (auto& p : transformed_points)
                {
                    tmp.nodes().add(p);
                }
                auto offset = polygon_element.connectivity().size();
                for (auto& s : boundary_segments)
                {
                    tmp.elements().add(ElementType::L1D2N, { s.first + offset, s.second + offset });
                }
                write_vtk(tmp, "failed_triangle_input.vtk", false);
            }

            auto global_p = simple_math::to_point(global_coord);
            total_mesh.nodes().add(global_p);
            nodes_in_triangulation.push_back(total_mesh.nodes().size() - 1);
        }
        // Now form the connectivities for the new triangles
        std::vector<LocalIndex> conn(3, 0);
        for (int i = 0; i < output.numberoftriangles; i++)
        {
            for (auto j : utility::range(3))
            {
                conn[j] = nodes_in_triangulation[output.trianglelist[i * 3 + j]];
            }
            total_mesh.elements().add(ElementType::T2D3N, conn);
            total_mesh.elements().back().component_id() = polygon_element.component_id();
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

    auto n_new_triangles = total_mesh.elements().size() - orig_n_elements;

    int adjust = 0;
    if (remove_polygon)
    {
        total_mesh.elements().remove_stable(polygon_element.local_index());
        adjust = 1;
    }

    std::vector<LocalIndex> new_triangle_indices;
    for (auto i : utility::range(n_new_triangles))
    {
        new_triangle_indices.push_back(orig_n_elements + i - adjust);
    }
    return new_triangle_indices;
}

std::pair<std::vector<LocalIndex>, std::vector<LocalIndex>> triangulate_all_polygons(
  Mesh& mesh,
  double max_tri_area)
{
    auto orig_n_elements = mesh.elements().size();

    // Loop over elements in surface mesh and submesh polygons with more than 4 nodes
    std::vector<LocalIndex> elements_to_remove;
    for (auto polygon_element_index : utility::range(mesh.elements().size()))
    {
        auto element = mesh.elements()[polygon_element_index];
        if (element.connectivity().size() > 3)
        {
            std::vector<std::vector<LocalIndex>> boundary_segments;
            auto conn = element.connectivity();
            for (auto segment : polygon_segments(element.connectivity()))
            {
                boundary_segments.push_back({ conn[segment.first], conn[segment.second] });
            }
            triangulate_polygon(mesh, element, max_tri_area, false);
            elements_to_remove.push_back(element.local_index());
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

    auto n_new_triangles = mesh.elements().size() - orig_n_elements;

    // Note: elements_to_remove is naturally sorted
    for (auto i : utility::range(elements_to_remove.size()))
    {
        mesh.elements().remove_stable(elements_to_remove[i] - static_cast<LocalIndex>(i));
    }


    std::vector<LocalIndex> new_triangle_indices;
    for (auto i : utility::range(n_new_triangles))
    {
        new_triangle_indices.push_back(orig_n_elements + i -
                                       static_cast<LocalIndex>(elements_to_remove.size()));
    }
    return { elements_to_remove, new_triangle_indices };
}
}
