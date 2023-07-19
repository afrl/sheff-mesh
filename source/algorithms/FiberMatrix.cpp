#include "mesh/algorithms/FiberMatrix.h"

#include "mesh/algorithms/Creation.h"
#include "mesh/algorithms/IO.h"
#include "mesh/algorithms/Polygons.h"
#include "simple_math/LineGeometry.h"
#include "simple_math/Utility.h"
#include "utility/Permutation.h"
#include "utility/EraseByIndex.h"

using namespace simple_math::point_operators;

namespace mesh
{
//////////////////////////////////////////////////////////////////////////////////////////////////
// Helper static functions

Mesh extract_boundary_mesh(const triangle::triangulateio& io)
{
    Mesh m;

    std::vector<bool> is_point_on_boundary(io.numberofpoints, false);
    for (auto i : utility::range(io.numberofedges))
    {
        if (io.edgemarkerlist[i] > 0)
        {
            is_point_on_boundary[io.edgelist[2 * i]]     = true;
            is_point_on_boundary[io.edgelist[2 * i + 1]] = true;
        }
    }

    std::vector<LocalIndex> input_to_new_point_map(io.numberofpoints, NOT_SET_LocalIndex);
    LocalIndex counter = 0;
    for (auto i : utility::range(io.numberofpoints))
    {
        if (is_point_on_boundary[i])
        {
            input_to_new_point_map[i] = counter++;
        }
    }

    for (auto i : utility::range(io.numberofpoints))
    {
        if (is_point_on_boundary[i])
        {
            m.nodes().add({ io.pointlist[2 * i], io.pointlist[2 * i + 1], 0. });
        }
    }
    for (auto i : utility::range(io.numberofedges))
    {
        if (io.edgemarkerlist[i] > 0)
        {
            m.elements().add(ElementType::L1D2N,
                             { input_to_new_point_map[io.edgelist[2 * i]],
                               input_to_new_point_map[io.edgelist[2 * i + 1]] });
            m.elements().back().component_id() = io.edgemarkerlist[i];
        }
    }
    return m;
}

Mesh extract_triangle_mesh(const triangle::triangulateio& io)
{
    Mesh m;
    for (auto i : utility::range(io.numberofpoints))
    {
        m.nodes().add({ io.pointlist[2 * i], io.pointlist[2 * i + 1], 0. });
    }
    for (auto i : utility::range(io.numberoftriangles))
    {
        m.elements().add(ElementType::T2D3N,
                         { (LocalIndex)io.trianglelist[3 * i],
                           (LocalIndex)io.trianglelist[3 * i + 1],
                           (LocalIndex)io.trianglelist[3 * i + 2] });
        if (io.numberoftriangleattributes > 0)
        {
            m.elements().back().component_id() = io.triangleattributelist[i];
        }
    }
    return m;
}

triangle::triangulateio triangulate(const Mesh& input,
                                    const std::vector<Point>& attribute_points,
                                    const std::vector<int>& attributes,
                                    const bool allow_new_nodes_on_boundary,
                                    const double max_triangle_area,
                                    double min_angle)
{
    _UTIL_ASSERT(attribute_points.size() == attributes.size(),
                 "Size mismatch betweeen attribute_points and attributes.");

    // Note - This method uses Triangle, a program by Jonathan Richard Shewchuk at Berkely
    // See his website at http://www.cs.cmu.edu/~quake/triangle.html for more information
    // This method is fashioned after tricall.c
    struct triangle::triangulateio in;

    in.numberofpoints          = static_cast<int>(input.nodes().size());
    in.numberofpointattributes = 0;
    in.pointlist               = new double[in.numberofpoints * 2];
    in.pointmarkerlist         = NULL;
    for (unsigned int i = 0; i < in.numberofpoints; ++i)
    {
        in.pointlist[i * 2]     = input.nodes()[i].x();
        in.pointlist[i * 2 + 1] = input.nodes()[i].y();
    }
    in.pointattributelist = NULL;
    in.numberofsegments   = static_cast<int>(input.elements().size());
    in.segmentlist        = new int[in.numberofsegments * 2];
    in.segmentmarkerlist  = new int[in.numberofsegments];
    for (int i = 0; i < in.numberofsegments; ++i)
    {
        in.segmentlist[i * 2]     = input.elements()[i].node(0).local_index();
        in.segmentlist[i * 2 + 1] = input.elements()[i].node(1).local_index();
        in.segmentmarkerlist[i]   = input.elements()[i].component_id();
    }
    in.numberofregions = attribute_points.size();
    in.regionlist      = new double[in.numberofregions * 4];
    for (unsigned int i = 0; i < attribute_points.size(); ++i)
    {
        in.regionlist[i * 4]     = attribute_points[i][0];
        in.regionlist[i * 4 + 1] = attribute_points[i][1];
        in.regionlist[i * 4 + 2] = attributes[i];
        in.regionlist[i * 4 + 3] = -1.0; // Area constraint, negative means none applied.
    }

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

    triangle::triangulateio output;
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

    triangulate(args, &in, &output, (struct triangle::triangulateio*)NULL);

    return output;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Definitions

Mesh make_rve_boundary(const std::vector<Point2D>& ordered_vertices, const double angle_tolerance)
{
    Mesh rve_boundary;

    for (auto i : utility::indices_as<LocalIndex>(ordered_vertices))
    {
        rve_boundary.nodes().add({ ordered_vertices[i][0], ordered_vertices[i][1], 0.0 });
    }

    auto segments = range_polygon_segments(ordered_vertices);
    for (auto& s : segments)
    {
        rve_boundary.elements().add(ElementType::L1D2N, utility::slice_c(s));
    }

    auto get_segment_normal = [&](const std::size_t i) -> simple_math::Vector2D
    {
        auto delta_x = ordered_vertices[segments[i][1]][0] - ordered_vertices[segments[i][0]][0];
        auto delta_y = ordered_vertices[segments[i][1]][1] - ordered_vertices[segments[i][0]][1];
        simple_math::Vector2D v(-delta_y, delta_x);
        v.normalize();
        return v;
    };

    // When the angle changes, mark it as a new component
    simple_math::Vector2D last_normal = get_segment_normal(0);
    ComponentIndex curr_marker        = 0;
    const double dot_tolerance        = cos(EIGEN_PI * angle_tolerance / 180.0);
    for (auto i : utility::indices(segments))
    {
        auto normal    = get_segment_normal(i);
        const auto dot = normal.dot(last_normal);
        if (dot < dot_tolerance)
        {
            curr_marker++;
        }
        last_normal = normal;

        rve_boundary.elements()[i].component_id() = curr_marker;
    }

    // Check that the last marker is not really just part of the first
    if (get_segment_normal(segments.size() - 1).dot(get_segment_normal(0)) >= dot_tolerance)
    {
        for (auto e : rve_boundary.elements())
        {
            if (e.component_id() == curr_marker)
            {
                e.component_id() = 0;
            }
        }
    }

    return rve_boundary;
}

struct RVEIntersection
{
    simple_math::CircleIntersection intersect;
    LocalIndex intersected_segment_index;
};

std::vector<RVEIntersection> intersect_fiber_rve_boundary(const StraightCircularFiber& fiber,
                                                          const Mesh& boundary)
{
    std::vector<RVEIntersection> intersections;
    for (auto& e : boundary.elements())
    {
        auto intersect = simple_math::intersect_circle_line_segment(fiber.radius,
                                                                    fiber.center,
                                                                    e.node(0).point(),
                                                                    e.node(1).point());
        if (intersect.type != simple_math::CircleIntersection::NO_INTERSECTION)
        {
            intersections.push_back({ intersect, e.local_index() });
        }
    }
    return intersections;
}

void sort_points_around_fiber(const StraightCircularFiber& fiber,
                              std::vector<Point>& points,
                              std::vector<bool>& is_point_on_boundary,
                              std::vector<LocalIndex>& intersected_rve_boundary_element)
{
    std::vector<std::pair<int, double>> quadrant_angle;
    for (auto& p : points)
    {
        auto offset  = p - fiber.center;
        int quadrant = 0;
        if (offset[0] >= 0.0 && offset[1] >= 0.0)
        {
            quadrant = 1;
        }
        else if (offset[0] < 0.0 && offset[1] >= 0.0)
        {
            quadrant = 2;
        }
        else if (offset[0] < 0.0 && offset[1] < 0.0)
        {
            quadrant = 3;
        }
        else
        {
            quadrant = 4;
        }

        double theta = 0.0;
        if (offset[0] == 0.0)
        {
            theta = simple_math::sign(offset[1]) * EIGEN_PI / 2.0;
        }
        else
        {
            theta = atan(offset[1] / offset[0]);
        }
        quadrant_angle.push_back({ quadrant, theta });
    }
    auto permuation = utility::find_sorted_permutation(
      quadrant_angle,
      [](const std::pair<int, double>& a, const std::pair<int, double>& b) -> bool
      {
          if (a.first == b.first)
          {
              return a.second < b.second;
          }
          else
          {
              return a.first < b.first;
          }
      });
    utility::apply_permutation(points, permuation);
    utility::apply_permutation(is_point_on_boundary, permuation);
    utility::apply_permutation(intersected_rve_boundary_element, permuation);
}

struct PointsAroundFiber
{
    std::vector<Point> points;
    std::vector<bool> is_point_on_boundary;
    std::vector<LocalIndex> intersected_rve_boundary_element;
};

PointsAroundFiber create_points_around_fiber(const StraightCircularFiber& fiber,
                                             const Mesh& rve_boundary,
                                             std::vector<Point> rve_boundary_points,
                                             const double target_edge_length)
{
    std::vector<Point> points;
    std::vector<bool> is_point_on_boundary;
    std::vector<LocalIndex> intersected_rve_boundary_element;

    // Find intersections with boundaries and include those points
    auto intersections = intersect_fiber_rve_boundary(fiber, rve_boundary);
    for (auto& i : intersections)
    {
        if (i.intersect.type != i.intersect.NO_INTERSECTION)
        {
            points.push_back(i.intersect.intersect_1);
            intersected_rve_boundary_element.push_back(i.intersected_segment_index);
        }
        if (i.intersect.type == i.intersect.TWO_INTERSECTS)
        {
            points.push_back(i.intersect.intersect_2);
            intersected_rve_boundary_element.push_back(i.intersected_segment_index);
        }
    }
    is_point_on_boundary.resize(points.size(), true);

    // Create evenly space points around the fiber
    const double fiber_perimeter       = 2.0 * EIGEN_PI * fiber.radius;
    int n_points_around_fiber          = fiber_perimeter / target_edge_length;
    double adjusted_target_edge_length = target_edge_length;
    if (n_points_around_fiber < 8)
    {
        n_points_around_fiber       = 8;
        adjusted_target_edge_length = fiber_perimeter / n_points_around_fiber;
    }

    auto n_intersection_points = points.size();
    for (auto i : utility::range(n_points_around_fiber))
    {
        double theta       = 2.0 * EIGEN_PI / n_points_around_fiber * i;
        const Point offset = { fiber.radius * cos(theta), fiber.radius * sin(theta), 0. };
        const Point point  = fiber.center + offset;
        if (simple_math::is_point_in_polygon(point, rve_boundary_points))
        {
            bool point_close_to_intersection = false;
            for (auto j : utility::range(n_intersection_points))
            {
                if (simple_math::distance(points[j], point) < 1e-2 * adjusted_target_edge_length)
                {
                    point_close_to_intersection = true;
                    break;
                }
            }
            if (!point_close_to_intersection)
            {
                points.push_back(point);
            }
        }
    }
    is_point_on_boundary.resize(points.size(), false);
    intersected_rve_boundary_element.resize(points.size(), 0);

    sort_points_around_fiber(fiber, points, is_point_on_boundary, intersected_rve_boundary_element);

    return { points, is_point_on_boundary, intersected_rve_boundary_element };
}

struct TriangleInput
{
    Mesh boundary_mesh;
    std::vector<Point> internal_fiber_points;
};

TriangleInput get_initial_triangle_input(const FiberMatrixRVEDescription& rve,
                                         const double target_edge_length)
{
    // Note: component IDs in triangle_mesher_input will be used to track the corresponding fiber or
    // RVE boundary segment. To do this, component IDs in the range [1, #fibers] correspond to
    // the fiber index + 1, and IDs in the range [#fibers + 1, #fibers + 1 + #RVE_boundary_segments]
    // correspond to the RVE boundary segments.  The output from Triangle will conveniently label
    // all other edges with an attribute of 0.

    Mesh triangle_mesher_input;
    for (auto& n : rve.rve_boundary.nodes())
    {
        triangle_mesher_input.nodes().add(n.point());
    }

    std::vector<Point> rve_boundary_points;
    for (auto& e : rve.rve_boundary.elements())
    {
        _UTIL_ASSERT(e.local_index() == e.node(0).local_index(), "");
        rve_boundary_points.push_back(e.node(0).point());
    }

    std::vector<Point> internal_fiber_point;
    std::vector<std::vector<LocalIndex>> boundary_segments_to_intersect_points_map(
      rve.rve_boundary.elements().size());
    for (auto fiber_i : utility::indices(rve.fibers))
    {
        auto& fiber = rve.fibers[fiber_i];

        auto result = create_points_around_fiber(fiber,
                                                 rve.rve_boundary,
                                                 rve_boundary_points,
                                                 target_edge_length);

        // Add nodes for intersection points
        // Note: it does track which intersection points lie on which RVE boundary segment
        auto node_offset = triangle_mesher_input.nodes().size();
        for (auto i : utility::indices(result.points))
        {
            triangle_mesher_input.nodes().add(result.points[i]);
            if (result.is_point_on_boundary[i])
            {
                const auto& new_point_index = triangle_mesher_input.nodes().back().local_index();
                const auto& intersected_element = result.intersected_rve_boundary_element[i];
                boundary_segments_to_intersect_points_map[intersected_element].push_back(
                  new_point_index);
            }
        }

        // Create segments around fiber (except boundaries)
        for (auto& s : utility::closed_range((LocalIndex)result.points.size()))
        {
            if (!result.is_point_on_boundary[s.i] || !result.is_point_on_boundary[s.j])
            {
                triangle_mesher_input.elements().add(ElementType::L1D2N,
                                                     { s.i + node_offset, s.j + node_offset });
                triangle_mesher_input.elements().back().component_id() = fiber_i + 1;
            }
        }

        internal_fiber_point.push_back(calc_mean_coordinate(result.points));
    }

    for (auto& e : rve.rve_boundary.elements())
    {
        // Note: copy intentional
        auto p_indices_for_segment = boundary_segments_to_intersect_points_map[e.local_index()];
        p_indices_for_segment.push_back(e.node(0).local_index());
        p_indices_for_segment.push_back(e.node(1).local_index());

        std::vector<Point> points_for_segment(p_indices_for_segment.size());
        for (auto i : utility::indices(points_for_segment))
        {
            points_for_segment[i] = triangle_mesher_input.nodes()[p_indices_for_segment[i]].point();
        }

        auto permutation = simple_math::find_sort_for_points_on_line(e.node(0).point(),
                                                                     e.node(1).point(),
                                                                     points_for_segment);
        utility::apply_permutation(p_indices_for_segment, permutation);

        for (auto& s : utility::open_range(permutation.size()))
        {
            triangle_mesher_input.elements().add(
              ElementType::L1D2N,
              { p_indices_for_segment[s.i], p_indices_for_segment[s.j] });
            triangle_mesher_input.elements().back().component_id() =
              e.component_id() + rve.fibers.size() + 1;
        }
    }

    return { triangle_mesher_input, internal_fiber_point };
}

void move_points_to_fibers(const FiberMatrixRVEDescription& rve, Mesh& boundary_mesh)
{
    // Note:
    //   fiber component IDs: [1, #fibers]
    //   RVE boundary component IDs: [#fibers + 1, #fibers + 1 + #RVE_boundary_segments]

    // Make sure not to move any points on the RVE boundary
    std::vector<bool> is_point_on_rve_boundary(boundary_mesh.nodes().size(), false);
    for (auto& e : boundary_mesh.elements())
    {
        if (e.component_id() > rve.fibers.size())
        {
            for (auto n : e.nodes())
            {
                is_point_on_rve_boundary[n.local_index()] = true;
            }
        }
    }

    // Adjust the coordinates of points that should be on the fiber
    for (auto& e : boundary_mesh.elements())
    {
        if (e.component_id() <= rve.fibers.size())
        {
            auto fiber_index = e.component_id() - 1;
            for (auto& n : e.nodes())
            {
                if (!is_point_on_rve_boundary[n.local_index()])
                {
                    Vector3D radial_vector =
                      n.coord() - simple_math::point_as_vector(rve.fibers[fiber_index].center);
                    radial_vector.normalize();

                    n.point() =
                      rve.fibers[fiber_index].center +
                      simple_math::to_point(rve.fibers[fiber_index].radius * radial_vector);
                }
            }
        }
    }
}

template<typename MeshType>
void get_refined_boundary(std::vector<std::vector<double>>& NodeCoords,
                          std::vector<std::vector<int>>& BoundarySegments,
                          std::vector<int>& SegmentMarkers,
                          const std::vector<std::vector<int>>& PeriodicMarkerPairs,
                          const std::vector<std::vector<double>>& PeriodicityVectors,
                          const int phase,
                          const double MaxTriangleArea,
                          double MinTriangleAngle)
{
    struct triangulateio mid;
    MeshCreation<MeshType>::triangle_kernel(mid,
                                            NodeCoords,
                                            BoundarySegments,
                                            SegmentMarkers,
                                            {},
                                            {},
                                            {},
                                            true,
                                            MaxTriangleArea,
                                            MinTriangleAngle);

    if (MeshCreationFiberMatrixDetail::output_debug_)
    {
        MeshType DebugMesh;
        auto TriAttributes = MeshCreation<MeshType>::get_mesh_from_triangle(DebugMesh, mid);
        MeshIO<MeshType>::write_VTK_mesh_file(DebugMesh,
                                              std::to_string(phase) + "_1_MeshToGetBoundary.vtk",
                                              false,
                                              false);
    }

    // Get the set of indices that lie on the boundaries
    std::set<int> refined_boundary_indices;
    for (auto i = 0; i < mid.numberofedges; i++)
        for (auto j = 0; j < 2; j++)
            refined_boundary_indices.insert(mid.edgelist[2 * i + j]);

    MeshType tmp_mesh;

    std::map<int, int> mid_to_refined_boundary_index_map;
    int new_index = 0;
    std::vector<NodePtr> refined_boundary_nodes;
    for (int index : refined_boundary_indices)
    {
        mid_to_refined_boundary_index_map[index] = new_index;
        std::vector<double> coord = { mid.pointlist[2 * index], mid.pointlist[2 * index + 1] };
        tmp_mesh.add_node_by_coord(coord);
        refined_boundary_nodes.push_back(tmp_mesh.nodes_back());
        new_index += 1;
    }

    // Build a map for getting the periodicity pair index given the marker
    std::map<int, size_t> marker_to_index_map;
    for (auto i = 0; i < PeriodicMarkerPairs.size(); i++)
    {
        marker_to_index_map[PeriodicMarkerPairs[i][0]] = i;
        marker_to_index_map[PeriodicMarkerPairs[i][1]] = i;
    }

    // Add all the segments for non-periodic and the first pair of the periodic boundaries
    std::vector<int> refined_edge_markers;
    for (auto i = 0; i < mid.numberofedges; i++)
    {
        int marker = mid.edgemarkerlist[i];
        if (marker != 0)
        {
            // Not a periodic marker
            if (marker_to_index_map.find(marker) == marker_to_index_map.end())
            {
                auto refined_index1 = mid_to_refined_boundary_index_map[mid.edgelist[2 * i]];
                auto refined_index2 = mid_to_refined_boundary_index_map[mid.edgelist[2 * i + 1]];
                std::vector<NodePtr> elem_nodes = { refined_boundary_nodes[refined_index1],
                                                    refined_boundary_nodes[refined_index2] };
                tmp_mesh.add_element();
                tmp_mesh.elements_back()->initialize(ElemGeomType::L1D2N, elem_nodes);
                refined_edge_markers.push_back(marker);
            }
            else
            {
                auto periodic_pair_index = marker_to_index_map[marker];
                // First of the pair will be added
                if (PeriodicMarkerPairs[periodic_pair_index][0] == marker)
                {
                    auto refined_index1 = mid_to_refined_boundary_index_map[mid.edgelist[2 * i]];
                    auto refined_index2 =
                      mid_to_refined_boundary_index_map[mid.edgelist[2 * i + 1]];
                    std::vector<NodePtr> elem_nodes = { refined_boundary_nodes[refined_index1],
                                                        refined_boundary_nodes[refined_index2] };
                    tmp_mesh.add_element();
                    tmp_mesh.elements_back()->initialize(ElemGeomType::L1D2N, elem_nodes);
                    refined_edge_markers.push_back(marker);

                    // Now add the second of the pair
                    // Add new nodes for the periodic segment, which will be removed later
                    std::vector<double> p1 = { mid.pointlist[2 * mid.edgelist[2 * i]],
                                               mid.pointlist[2 * mid.edgelist[2 * i] + 1] };
                    auto n1                = tmp_mesh.add_node_by_coord(
                      add_vec(p1, PeriodicityVectors[periodic_pair_index]));
                    std::vector<double> p2 = { mid.pointlist[2 * mid.edgelist[2 * i + 1]],
                                               mid.pointlist[2 * mid.edgelist[2 * i + 1] + 1] };
                    auto n2                = tmp_mesh.add_node_by_coord(
                      add_vec(p2, PeriodicityVectors[periodic_pair_index]));
                    // Note: reverse connectivity order since it is the other side
                    std::vector<NodePtr> elem_nodes2 = { n2, n1 };
                    tmp_mesh.add_element();
                    tmp_mesh.elements_back()->initialize(ElemGeomType::L1D2N, elem_nodes2);
                    refined_edge_markers.push_back(PeriodicMarkerPairs[periodic_pair_index][1]);
                }
            }
        }
    }

    // frees memory used by triangle
    free(mid.pointlist);
    free(mid.pointattributelist);
    free(mid.pointmarkerlist);
    free(mid.trianglelist);
    free(mid.triangleattributelist);
    free(mid.neighborlist);
    free(mid.segmentlist);
    free(mid.segmentmarkerlist);
    free(mid.edgelist);
    free(mid.edgemarkerlist);

    // Now remove unattached and duplicate nodes
    MeshManipulation<MeshType>::remove_unattached_nodes(tmp_mesh, 0);
    MeshManipulation<MeshType>::remove_duplicate_nodes(tmp_mesh);
    tmp_mesh.update_global_numbering();

    if (MeshCreationFiberMatrixDetail::output_debug_)
        MeshIO<MeshType>::write_VTK_mesh_file(tmp_mesh,
                                              std::to_string(phase) +
                                                "_2_ExtractedBoundaryMesh.vtk",
                                              false,
                                              false);

    // Now extract the boundary info
    NodeCoords.clear();
    for (auto n_iter = tmp_mesh.nodes_begin(); n_iter != tmp_mesh.nodes_end(); n_iter++)
        NodeCoords.push_back((*n_iter)->get_coord());
    BoundarySegments.clear();
    for (auto e_iter = tmp_mesh.elements_begin(); e_iter != tmp_mesh.elements_end(); e_iter++)
    {
        std::vector<int> segment;
        for (auto ind : (*e_iter)->get_local_conn())
            segment.push_back(ind);
        BoundarySegments.push_back(segment);
    }

    SegmentMarkers = refined_edge_markers;
}

Mesh mesh_rve(const FiberMatrixRVEDescription& rve, const double target_edge_length)
{
    // Note: component IDs in triangle_mesher_input will be used to track the corresponding fiber or
    // RVE boundary segment. To do this, component IDs in the range [1, #fibers] correspond to
    // the fiber index + 1, and IDs in the range [#fibers + 1, #fibers + 1 + #RVE_boundary_segments]
    // correspond to the RVE boundary segments.  The output from Triangle will conveniently label
    // all other edges with an attribute of 0.
    auto triangle_mesher_input = get_initial_triangle_input(rve, target_edge_length);

    write_vtk(triangle_mesher_input.boundary_mesh, "triangle_mesher_input.vtk");

    std::vector<int> fiber_ids(rve.fibers.size());
    std::iota(fiber_ids.begin(), fiber_ids.end(), 1);
    auto triangulation = triangulate(triangle_mesher_input.boundary_mesh,
                                     triangle_mesher_input.internal_fiber_points,
                                     fiber_ids,
                                     true,
                                     0.5 * target_edge_length * target_edge_length,
                                     25.0);

    write_vtk(extract_triangle_mesh(triangulation), "output_triangles.vtk");

    auto boundary_mesh = extract_boundary_mesh(triangulation);
    write_vtk(boundary_mesh, "output_boundary.vtk");
    move_points_to_fibers(rve, boundary_mesh);
    write_vtk(boundary_mesh, "adjusted_boundary.vtk");

    // TODO transfer points across periodic boundaries
    //
    // Note:
    //   fiber component IDs: [1, #fibers]
    //   RVE boundary component IDs: [#fibers + 1, #fibers + 1 + #RVE_boundary_segments]

    // We will need to define a line to sort the points along
    std::vector<std::array<Point, 2>> rve_boundary_marker_to_line_points_map(
      rve.rve_boundary.n_components());


    std::vector<std::set<LocalIndex>> rve_boundary_marker_to_points_map(
      rve.rve_boundary.n_components());
    std::vector<std::vector<LocalIndex>> rve_boundary_marker_to_elements_map(
      rve.rve_boundary.n_components());
    auto boundary_component_id_start = rve.fibers.size() + 1 + rve.rve_boundary.n_components();
    for (auto& e : boundary_mesh.elements())
    {
        if (e.component_id() >= boundary_component_id_start)
        {
            const auto rve_boundary_marker = e.component_id() - boundary_component_id_start;
            for (auto& n : e.nodes())
            {
                rve_boundary_marker_to_points_map[rve_boundary_marker].insert(n.local_index());
            }
            rve_boundary_marker_to_elements_map[rve_boundary_marker].push_back(e.local_index());
        }
    }

    std::vector<bool> is_point_on_rve_boundary(boundary_mesh.nodes().size(), false);
    std::vector<bool> is_point_on_fiber(boundary_mesh.nodes().size(), false);
    for (auto& e : boundary_mesh.elements())
    {
        if (e.component_id() > rve.fibers.size())
        {
            for (auto n : e.nodes())
            {
                is_point_on_rve_boundary[n.local_index()] = true;
            }
        }
        else if (e.component_id() >= 1 && e.component_id() <= rve.fibers.size())
        {
            for (auto n : e.nodes())
            {
                is_point_on_fiber[n.local_index()] = true;
            }
        }
    }

    for (auto i : utility::indices(rve.periodic_markers))
    {
        auto& periodic_marker_pair = rve.periodic_markers[i];
        auto& periodicity_vector   = rve.periodicity_vectors[i];

        // Create a list of points that includes those on the first marker and a projection of those
        // on the second marker.
        std::vector<Point> points;
        std::vector<bool> is_point_on_intersection;
        for (auto& point_index : rve_boundary_marker_to_points_map[periodic_marker_pair[0]])
        {
            points.push_back(boundary_mesh.nodes()[point_index].point());
            is_point_on_intersection.push_back(is_point_on_fiber[point_index] &&
                                               is_point_on_rve_boundary[point_index]);
        }
        for (auto& point_index : rve_boundary_marker_to_points_map[periodic_marker_pair[1]])
        {
            points.push_back(boundary_mesh.nodes()[point_index].point() - periodicity_vector);
            is_point_on_intersection.push_back(is_point_on_fiber[point_index] &&
                                               is_point_on_rve_boundary[point_index]);
        }

        // Sort the points along the line
        auto element_on_marker_1 = rve_boundary_marker_to_elements_map[periodic_marker_pair[0]][0];
        auto p1_on_marker_1      = boundary_mesh.elements()[element_on_marker_1].node(0).point();
        auto p2_on_marker_1      = boundary_mesh.elements()[element_on_marker_1].node(1).point();
        auto p = simple_math::find_sort_for_points_on_line(p1_on_marker_1, p2_on_marker_1, points);
        utility::apply_permutation(points, p);
        utility::apply_permutation(is_point_on_intersection, p);

        // Remove any points within a distance tolerance to the previous point to avoid duplicates
        std::vector<std::size_t> points_to_erase;
        std::size_t last_non_erased_j = 0;
        for (auto j : utility::range_as<std::size_t>(1, points.size() - 1))
        {
            if (!is_point_on_intersection[j])
            {
                if (simple_math::distance(points[j], points[last_non_erased_j]) <
                    1e-4 * target_edge_length)
                {
                    points_to_erase.push_back(j);
                }
                else
                {
                    last_non_erased_j = j;
                }
            }
            else
            {
                last_non_erased_j = j;
            }
        }
        utility::erase_by_indices(points, points_to_erase);
        utility::erase_by_indices(is_point_on_intersection, points_to_erase);
        points_to_erase.clear();

        last_non_erased_j = points.size() - 1;
        for (auto j : utility::reversed_range_as<std::size_t>(1, points.size() - 1))
        {
            if (!is_point_on_intersection[j])
            {
                if (simple_math::distance(points[j], points[last_non_erased_j]) <
                    1e-4 * target_edge_length)
                {
                    points_to_erase.push_back(j);
                }
                else
                {
                    last_non_erased_j = j;
                }
            }
            else
            {
                last_non_erased_j = j;
            }
        }
        utility::erase_by_indices(points, points_to_erase);
        utility::erase_by_indices(is_point_on_intersection, points_to_erase);
        points_to_erase.clear();
        // Remove points and Sweep opposite way
    }

    // TODO track which nodes/segments lie on which fibers/boundary

    Mesh m;
    return m;
}
}