#include "mesh/algorithms/Polygons.h"

#include "second_party_mesh/triangle/triangle.h"
#include "simple_math/Utility.h"

using namespace simple_math::point_operators;

namespace mesh
{
//////////////////////////////////////////////////////////////////////////////////////////////////
// Helpers

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

SegmentStats calc_segment_length_stats(const Mesh& m)
{
    std::vector<double> segment_lengths;
    for (auto& element : m.elements())
    {
        auto& conn = element.connectivity();
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
    auto& conn = polygon_element.connectivity();
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

void triangulate_polygon(Mesh& total_mesh,
                         const ElementConstAccessor& polygon_element,
                         double max_area)
{
    assert(&polygon_element.node_container() == &total_mesh.nodes());

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
        }
    }

    total_mesh.elements().remove_stable(polygon_element.local_index());

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
}