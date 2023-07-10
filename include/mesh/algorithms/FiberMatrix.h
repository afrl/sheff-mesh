#pragma once

// BetaMesh dependencies
// #include "algorithms/GroupOperations.h"
// #include "algorithms/MeshCreation.h"
// #include "algorithms/MeshManipulation.h"
// #include "algorithms/MeshSearch.h"
// #include "misc/BetaMeshUtility.h"

#include "mesh/Mesh.h"

#include <iostream>
#include <set>
#include <vector>

namespace mesh
{

struct StraightCircularFiber
{
    Point center;
    double radius;
};

struct FiberMatrixRVEDescription
{
    /// List of fibers in the RVE (note: right now it is just circular straight fibers, but this
    /// could be extended)
    std::vector<StraightCircularFiber> fibers;

    /// A mesh describing the boundary of the entire RVE discretized as a set of line segments that
    /// form a closed polygon.
    Mesh rve_boundary;

    /// A list of all pairs of markers (component IDs) that should be periodic in the discretization
    std::vector<std::array<ComponentIndex, 2>> periodic_markers;
    /// The vector going from the first to the second marker for each pair of periodic markers.
    std::vector<Point> periodicity_vectors;
};

/// @brief A helper function to create the discretization for an RVE boundary based on the list of
///        ordered points.
///
/// @param angle_tolerance - starts a new facet if the angle whenever the angle between two adjacent
///        segments is larger than this tolerance.
Mesh make_rve_boundary(const std::vector<Point2D>& ordered_vertices, const double angle_tolerance);

Mesh mesh_rve(const FiberMatrixRVEDescription& rve, const double target_edge_length);



/*

/// A namespace for static variables that are shared by all templated
/// versions of MeshCreationFiberMatrix.
class MeshCreationFiberMatrixDetail
{
  public:
    /// A bool indicating if intermediate meshes are output during the meshing process.
    static bool output_debug_;
};

/// This collection of mesh algorithms pertains to creating meshes for randomly generated
/// fiber/matrix RVE's.
/// \attention Template requirement: MeshType must be Mesh type or a derived class.
template<typename MeshType>
class MeshCreationFiberMatrix : public MeshCreation<MeshType>
{
  public:
    typedef typename MeshType::NodeType NodeType;
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::ElementType ElementType;
    typedef typename MeshType::NodePtr NodePtr;
    typedef typename MeshType::FacePtr FacePtr;
    typedef typename MeshType::ElementPtr ElementPtr;

    /// Create a fiber matrix mesh
    /// \param NodeCoords - the list of nodal coordinates
    /// \param BoundarySegments - the list of boundary segments that defines the RVE boundaries
    ///   and fiber/matrix boundaries (the integers refer to indices in the NodeCoords parameter)
    /// \param SegmentMarkers - the list of markers for all the boundaries.  Positve markers
    /// indicate
    ///   a fiber/matrix boundary, where each fiber has a unique value and starts at 0.  Negative
    ///   markers indciate an RVE boundary, where each boundary has a unique value and starts at 0.
    /// \param PeriodicMarkerPairs - the list of marker pairs (see SegmentMarkers) that are
    /// periodic.
    ///   each vector must be a length 2. (Vector is used for Python compatibility)
    /// \param PeriodicityVectors - vector of periodicity for each entry in PeriodicMarkerPairs.
    ///   The implicit order mastches PeriodicMarkerPairs.
    /// \param FiberCenters - the centers of all the fibers (the center may lie outside the RVE
    /// bounds) \param FiberRadii - the radii of all of the fibers (the length must match
    /// FiberCenters) \param FiberInternalCoords - a coordinate that lies in the fiber and inside
    /// the RVE, and like
    ///   FiberRadii, the lenght must match FiberCenters.
    /// \param MaxTriangleArea - the maximum area a triangle is allowed to have (may be violated for
    /// some
    ///   triangles lying on the boundary of the RVE and fiber)
    /// \param MinTriangleAngle - the minimum angle between any two triangle edges (may be violated
    /// for some
    ///   triangles lying on the boundary of the RVE and fiber)
    static Mesh make_fiber_matrix_rve_mesh(
      const std::vector<std::vector<double>>& NodeCoords,
      const std::vector<std::vector<int>>& BoundarySegments,
      const std::vector<int>& SegmentMarkers,
      const std::vector<std::vector<int>>& PeriodicMarkerPairs,
      const std::vector<std::vector<double>>& PeriodicityVectors,
      const std::vector<std::vector<double>>& FiberCenters,
      const std::vector<double>& FiberRadii,
      const std::vector<std::vector<double>>& FiberInternalCoords,
      const double MaxTriangleArea = -1.0,
      double MinTriangleAngle      = 28.6);

    /// A utility method for fiber/matrix RVE meshes.
    /// Returns a vertices object needed for plotting the triangles as polygons
    /// in matplotlib
    static std::vector<std::vector<std::vector<double>>> get_plot_polys(
      MeshType& mesh,
      const unsigned int MaxNodes = 100000);

    /// A utility method for fiber/matrix RVE meshes.
    /// Returns a color code object needed for plotting the triangles as polygons
    /// in matplotlib with the color corresponding to domain id
    static std::vector<std::vector<double>> get_plot_poly_colors(
      MeshType& mesh,
      const unsigned int MaxNodes = 100000);

    /// The FiberMatrixModelGenerator tool will generate a .mesh_input file, which is read
    /// by this method and meshed.  It is useful for debugging or loading a past RVE
    static void read_fm_mesher_input_file(MeshType& mesh, std::string mesherInputFilepath);

    /// Sets the flag indicating if intermediate meshes are output during the meshing process.
    /// Note: this flag applies to all templated versions of MeshCreationFiberMatrix.
    static void set_output_debug_flag(bool flag)
    {
        MeshCreationFiberMatrixDetail::output_debug_ = flag;
    };

  protected:
    /// Given a coarse definition of the RVE discretization, this method
    /// fills the data with a refined discretization based on Triangle's
    /// refinement algorithms.
    static void get_refined_boundary(std::vector<std::vector<double>>& NodeCoords,
                                     std::vector<std::vector<int>>& BoundarySegments,
                                     std::vector<int>& SegmentMarkers,
                                     const std::vector<std::vector<int>>& PeriodicMarkerPairs,
                                     const std::vector<std::vector<double>>& PeriodicityVectors,
                                     const int phase,
                                     const double MaxTriangleArea = -1.0,
                                     double MinTriangleAngle      = 28.6);

    /// Adjust nodes that lie on the fiber/matrix boundaries to lie right
    /// on the fiber/matrix boundary.
    static void adjust_fiber_nodes(std::vector<std::vector<double>>& NodeCoords,
                                   const std::vector<std::vector<int>>& BoundarySegments,
                                   const std::vector<int>& SegmentMarkers,
                                   const std::vector<std::vector<double>>& FiberCenters,
                                   const std::vector<double>& FiberRadii);

    /// Given a good boundary definition, this creates the final quadratic mesh.
    static void create_refined_quadratic_mesh(
      MeshType& mesh,
      std::vector<std::vector<double>>& NodeCoords,
      std::vector<std::vector<int>>& BoundarySegments,
      std::vector<int>& SegmentMarkers,
      const std::vector<std::vector<double>>& FiberCenters,
      const std::vector<double>& FiberRadii,
      const std::vector<std::vector<double>>& FiberInternalCoords,
      const double MaxTriangleArea = -1.0,
      double MinTriangleAngle      = 28.6);

    /// Given the nodal coordinates and boundary segments, this creates a mesh
    /// of 2 node line elements for visualization purposes.
    static MeshType extract_boundary_mesh(std::vector<std::vector<double>>& NodeCoords,
                                          std::vector<std::vector<int>>& BoundarySegments,
                                          std::vector<int>& SegmentMarkers);
};

template<typename MeshType>
inline void MeshCreationFiberMatrix<MeshType>::make_fiber_matrix_rve_mesh(
  MeshType& mesh,
  const std::vector<std::vector<double>>& NodeCoords,
  const std::vector<std::vector<int>>& BoundarySegments,
  const std::vector<int>& SegmentMarkers,
  const std::vector<std::vector<int>>& PeriodicMarkerPairs,
  const std::vector<std::vector<double>>& PeriodicityVectors,
  const std::vector<std::vector<double>>& FiberCenters,
  const std::vector<double>& FiberRadii,
  const std::vector<std::vector<double>>& FiberInternalCoords,
  const double MaxTriangleArea,
  double MinTriangleAngle)
{
    auto node_coords_copy       = NodeCoords;
    auto boundary_segments_copy = BoundarySegments;
    auto segment_markers_copy   = SegmentMarkers;

    if (MeshCreationFiberMatrixDetail::output_debug_)
    {
        MeshType tmp_mesh =
          MeshCreationFiberMatrix<MeshType>::extract_boundary_mesh(node_coords_copy,
                                                                   boundary_segments_copy,
                                                                   segment_markers_copy);
        MeshIO<MeshType>::write_VTK_mesh_file(tmp_mesh, "0_ExternalInput.vtk", false, false);
    }

    double min_radius               = *std::min_element(FiberRadii.begin(), FiberRadii.end());
    auto max_triangle_area_adjusted = M_PI * min_radius * min_radius / 20.0;
    // std::cout << "minarea " << min_radius << " " << max_triangle_area_adjusted << "\n";
    if (MaxTriangleArea < max_triangle_area_adjusted && MaxTriangleArea > 0.0)
        max_triangle_area_adjusted = MaxTriangleArea;

    // Let Triangle create nodes along the boundaries, extract them, duplicate for
    // periodicity, and adjust nodes to lie on the fibers
    {
        get_refined_boundary(node_coords_copy,
                             boundary_segments_copy,
                             segment_markers_copy,
                             PeriodicMarkerPairs,
                             PeriodicityVectors,
                             1,
                             max_triangle_area_adjusted,
                             MinTriangleAngle);

        adjust_fiber_nodes(node_coords_copy,
                           boundary_segments_copy,
                           segment_markers_copy,
                           FiberCenters,
                           FiberRadii);
    }

    // Let Triangle create nodes along the boundaries one more time since the periodicity
    // vectors can create very bad triangles unless new nodes allowed to be created on the
    // fibers lying on periodic boundaries
    {
        get_refined_boundary(node_coords_copy,
                             boundary_segments_copy,
                             segment_markers_copy,
                             PeriodicMarkerPairs,
                             PeriodicityVectors,
                             2,
                             max_triangle_area_adjusted,
                             MinTriangleAngle);

        adjust_fiber_nodes(node_coords_copy,
                           boundary_segments_copy,
                           segment_markers_copy,
                           FiberCenters,
                           FiberRadii);
    }

    // For degbugging, create a mesh for boundaries.
    if (MeshCreationFiberMatrixDetail::output_debug_)
    {
        MeshType tmp_mesh =
          MeshCreationFiberMatrix<MeshType>::extract_boundary_mesh(node_coords_copy,
                                                                   boundary_segments_copy,
                                                                   segment_markers_copy);
        MeshIO<MeshType>::write_VTK_mesh_file(tmp_mesh, "3_AdjustedBoundaries.vtk", false, false);
    }

    create_refined_quadratic_mesh(mesh,
                                  node_coords_copy,
                                  boundary_segments_copy,
                                  segment_markers_copy,
                                  FiberCenters,
                                  FiberRadii,
                                  FiberInternalCoords,
                                  max_triangle_area_adjusted,
                                  MinTriangleAngle);
}

template<typename MeshType>
inline std::vector<std::vector<std::vector<double>>>
MeshCreationFiberMatrix<MeshType>::get_plot_polys(MeshType& mesh, const unsigned int MaxNodes)
{
    std::vector<std::vector<std::vector<double>>> vertices;
    if (mesh.nodes_size() < MaxNodes)
    {
        for (auto e_iter = mesh.elements_begin(); e_iter != mesh.elements_end(); e_iter++)
        {
            std::vector<std::vector<double>> triangle_verts;
            for (auto n : (*e_iter)->get_conn())
                triangle_verts.push_back(n->get_coord());
            vertices.push_back(triangle_verts);
        }
    }
    return vertices;
}

template<typename MeshType>
inline std::vector<std::vector<double>> MeshCreationFiberMatrix<MeshType>::get_plot_poly_colors(
  MeshType& mesh,
  const unsigned int MaxNodes)
{
    std::vector<std::vector<double>> colors;
    if (mesh.nodes_size() < MaxNodes)
    {
        for (auto e_iter = mesh.elements_begin(); e_iter != mesh.elements_end(); e_iter++)
        {
            auto domainID = (*e_iter)->get_domain_id();
            if (domainID == 1) // Fibers
                colors.push_back({ 0.0, 0.0, 1.0, 1.0 });
            else
                colors.push_back({ 1.0, 0.0, 0.0, 0.5 });
        }
    }
    return colors;
}

template<typename MeshType>
inline void MeshCreationFiberMatrix<MeshType>::read_fm_mesher_input_file(
  MeshType& mesh,
  std::string mesherInputFilepath)
{
    ifstream myfile;
    myfile.open(mesherInputFilepath);
    int numPoints;
    myfile >> numPoints;
    std::vector<std::vector<double>> points;
    for (int i = 0; i < numPoints; i++)
    {
        std::vector<double> point(2, 0.0);
        myfile >> point[0] >> point[1];
        points.push_back(point);
    }
    int numSegments;
    myfile >> numSegments;
    std::vector<std::vector<int>> segments;
    for (int i = 0; i < numSegments; i++)
    {
        std::vector<int> segment(2, 0);
        myfile >> segment[0] >> segment[1];
        segments.push_back(segment);
    }
    int numSegmentMarkers;
    myfile >> numSegmentMarkers;
    std::vector<int> segmentMarkers(numSegmentMarkers, 0);
    for (int i = 0; i < numSegmentMarkers; i++)
    {
        myfile >> segmentMarkers[i];
    }
    int numPeriodicPairs;
    myfile >> numPeriodicPairs;
    std::vector<std::vector<int>> periodicPairs;
    for (int i = 0; i < numPeriodicPairs; i++)
    {
        int n1, n2;
        myfile >> n1 >> n2;
        periodicPairs.push_back({ n1, n2 });
    }
    std::vector<std::vector<double>> periodicVectors;
    for (int i = 0; i < numPeriodicPairs; i++)
    {
        double n1, n2;
        myfile >> n1 >> n2;
        periodicVectors.push_back({ n1, n2 });
    }
    int numFibers;
    myfile >> numFibers;
    std::vector<std::vector<double>> fiberCenters;
    for (int i = 0; i < numFibers; i++)
    {
        double n1, n2;
        myfile >> n1 >> n2;
        fiberCenters.push_back({ n1, n2 });
    }
    std::vector<double> fiberRadii;
    for (int i = 0; i < numFibers; i++)
    {
        double n1;
        myfile >> n1;
        fiberRadii.push_back(n1);
    }
    std::vector<std::vector<double>> fiberInternalCoords;
    for (int i = 0; i < numFibers; i++)
    {
        double n1, n2;
        myfile >> n1 >> n2;
        fiberInternalCoords.push_back({ n1, n2 });
    }
    myfile.close();
    MeshCreationFiberMatrix<MeshType>::make_fiber_matrix_rve_mesh(mesh,
                                                                  points,
                                                                  segments,
                                                                  segmentMarkers,
                                                                  periodicPairs,
                                                                  periodicVectors,
                                                                  fiberCenters,
                                                                  fiberRadii,
                                                                  fiberInternalCoords);
}

template<typename MeshType>
inline void MeshCreationFiberMatrix<MeshType>::get_refined_boundary(
  std::vector<std::vector<double>>& NodeCoords,
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

template<typename MeshType>
inline void MeshCreationFiberMatrix<MeshType>::adjust_fiber_nodes(
  std::vector<std::vector<double>>& NodeCoords,
  const std::vector<std::vector<int>>& BoundarySegments,
  const std::vector<int>& SegmentMarkers,
  const std::vector<std::vector<double>>& FiberCenters,
  const std::vector<double>& FiberRadii)
{
    for (auto seg_ind = 0; seg_ind < SegmentMarkers.size(); seg_ind++)
    {
        // If this is a fiber boundary segment
        if (SegmentMarkers[seg_ind] > 0)
        {
            auto fiber_ind    = SegmentMarkers[seg_ind] - 1; // Markers were 1 based
            auto fiber_center = FiberCenters[fiber_ind];
            auto fiber_radius = FiberRadii[fiber_ind];
            // Adjust both nodes of the segment (the nodes will be adjusted multiple times,
            // but it is fine.
            for (auto coord_ind : BoundarySegments[seg_ind])
            {
                auto dir_vec = subtract_vec(NodeCoords[coord_ind], fiber_center);
                normalize(dir_vec);
                NodeCoords[coord_ind] = add_vec(fiber_center, scale_vec(dir_vec, fiber_radius));
            }
        }
    }
}

template<typename MeshType>
inline void MeshCreationFiberMatrix<MeshType>::create_refined_quadratic_mesh(
  MeshType& mesh,
  std::vector<std::vector<double>>& NodeCoords,
  std::vector<std::vector<int>>& BoundarySegments,
  std::vector<int>& SegmentMarkers,
  const std::vector<std::vector<double>>& FiberCenters,
  const std::vector<double>& FiberRadii,
  const std::vector<std::vector<double>>& FiberInternalCoords,
  const double MaxTriangleArea,
  double MinTriangleAngle)
{
    // Note: this algorithm has been changed because of an important detail of the Triangle mesher
    // Triangle adds new nodes to the fiber/matrix boundaries because they are "internal" despite
    // the the flag to not add new nodes.  Consequently, the meshes for the matrix and the fibers
    // (collectively) are generated separately and then combined at the end.

    if (MeshCreationFiberMatrixDetail::output_debug_)
    {
        MeshType BoundaryMesh =
          MeshCreationFiberMatrix<MeshType>::extract_boundary_mesh(NodeCoords,
                                                                   BoundarySegments,
                                                                   SegmentMarkers);
        MeshIO<MeshType>::write_VTK_mesh_file(BoundaryMesh, "4_TotalBoundary.vtk", false, false);

        std::ofstream internalCoordFile("4_2_InternalCoords.txt");
        internalCoordFile << "x_coord"
                          << ","
                          << "y_coord" << std::endl;
        for (auto& coord : FiberInternalCoords)
            internalCoordFile << coord[0] << ", " << coord[1] << std::endl;
        internalCoordFile.close();
    }

    // Create the refined linear element mesh for the matrix
    {
        struct triangulateio mid;
        MeshCreation<MeshType>::triangle_kernel(mid,
                                                NodeCoords,
                                                BoundarySegments,
                                                SegmentMarkers,
                                                FiberInternalCoords,
                                                {},
                                                {},
                                                false,
                                                MaxTriangleArea,
                                                MinTriangleAngle);

        auto TriAttributes = MeshCreation<MeshType>::get_mesh_from_triangle(mesh, mid);
        MeshManipulation<MeshType>::make_mesh_quadratic(mesh);
        mesh.update_global_numbering();
        if (MeshCreationFiberMatrixDetail::output_debug_)
            MeshIO<MeshType>::write_VTK_mesh_file(mesh, "5_MatrixMesh.vtk", false, false);

        // Adjust midside nodes
        std::map<size_t, int> edge_hash_to_fiber_index;
        for (auto edge_ind = 0; edge_ind < mid.numberofedges; edge_ind++)
        {
            // If this is a fiber boundary segment
            if (mid.edgemarkerlist[edge_ind] > 0)
            {
                std::vector<int> fiber_boundary_edge = { mid.edgelist[2 * edge_ind],
                                                         mid.edgelist[2 * edge_ind + 1] };
                std::sort(fiber_boundary_edge.begin(), fiber_boundary_edge.end());
                auto hashID = VectorValueHasher<int>::hash(fiber_boundary_edge);
                // Note: the index is 0 based and the marker is 1 based
                edge_hash_to_fiber_index[hashID] = mid.edgemarkerlist[edge_ind] - 1;
            }
        }
        for (auto e_iter = mesh.elements_begin(); e_iter != mesh.elements_end(); e_iter++)
        {
            for (auto& f : (*e_iter)->get_faces())
            {
                int num1                                 = f->get_conn()[0]->get_global_number();
                int num2                                 = f->get_conn()[2]->get_global_number();
                std::vector<int> edge_linearized_indices = { num1, num2 };
                std::sort(edge_linearized_indices.begin(), edge_linearized_indices.end());
                auto hashID    = VectorValueHasher<int>::hash(edge_linearized_indices);
                auto find_iter = edge_hash_to_fiber_index.find(hashID);
                if (find_iter != edge_hash_to_fiber_index.end())
                {
                    auto fiber_ind    = find_iter->second;
                    auto fiber_center = FiberCenters[fiber_ind];
                    auto fiber_radius = FiberRadii[fiber_ind];
                    auto& face_conn   = f->get_conn();
                    for (auto i = 0; i < face_conn.size(); i++)
                    {
                        auto coord   = face_conn[i]->get_coord();
                        auto dir_vec = subtract_vec(coord, fiber_center);
                        normalize(dir_vec);
                        coord           = add_vec(fiber_center, scale_vec(dir_vec, fiber_radius));
                        (*face_conn[i]) = coord;
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
    }

    if (mesh.elements_size() == 0)
        throw _BETA_MSG_EXCEPT2(
          "Something caused the matrix mesh fail in triangle.  If this continues to happen, call "
          "MeshCreationFiberMatrix::set_output_debug_flag(bool) before creating the mesh to get "
          "more information.");

    // Create mesh for fibers
    {
        // If a fiber happens to straddle two boundaries, such as:
        //  | __
        //  |/  \
            //  |    |
        //  L\__/___
        //  ^ (this is the problem corner
        // There will not be a hole where needed to get rid of the matrix elements
        // A better way to do this is to just add boundary segments that define all the fibers, but
        // fow now, we will just make a whole lot of holes
        std::vector<std::vector<double>> matrixHoles;
        for (auto& element : mesh.elements())
            matrixHoles.push_back(element->get_internal_coordinate());

        MeshType FibersMesh;
        struct triangulateio mid;
        MeshCreation<MeshType>::triangle_kernel(mid,
                                                NodeCoords,
                                                BoundarySegments,
                                                SegmentMarkers,
                                                matrixHoles,
                                                {},
                                                {},
                                                false,
                                                MaxTriangleArea,
                                                MinTriangleAngle);

        auto TriAttributes = MeshCreation<MeshType>::get_mesh_from_triangle(FibersMesh, mid);
        MeshManipulation<MeshType>::make_mesh_quadratic(FibersMesh);

        FibersMesh.update_global_numbering();
        // Set domain IDs for fibers to 1
        GroupOperations<MeshType>::set_all_domain_ids(FibersMesh, 1);
        if (MeshCreationFiberMatrixDetail::output_debug_)
            MeshIO<MeshType>::write_VTK_mesh_file(FibersMesh, "6_FiberMesh.vtk", false, false);

        // Adjust midside nodes
        std::map<size_t, int> edge_hash_to_fiber_index;
        for (auto edge_ind = 0; edge_ind < mid.numberofedges; edge_ind++)
        {
            // If this is a fiber boundary segment
            if (mid.edgemarkerlist[edge_ind] > 0)
            {
                std::vector<int> fiber_boundary_edge = { mid.edgelist[2 * edge_ind],
                                                         mid.edgelist[2 * edge_ind + 1] };
                std::sort(fiber_boundary_edge.begin(), fiber_boundary_edge.end());
                auto hashID = VectorValueHasher<int>::hash(fiber_boundary_edge);
                // Note: the index is 0 based and the marker is 1 based
                edge_hash_to_fiber_index[hashID] = mid.edgemarkerlist[edge_ind] - 1;
            }
        }
        for (auto e_iter = FibersMesh.elements_begin(); e_iter != FibersMesh.elements_end();
             e_iter++)
        {
            for (auto& f : (*e_iter)->get_faces())
            {
                int num1                                 = f->get_conn()[0]->get_global_number();
                int num2                                 = f->get_conn()[2]->get_global_number();
                std::vector<int> edge_linearized_indices = { num1, num2 };
                std::sort(edge_linearized_indices.begin(), edge_linearized_indices.end());
                auto hashID    = VectorValueHasher<int>::hash(edge_linearized_indices);
                auto find_iter = edge_hash_to_fiber_index.find(hashID);
                if (find_iter != edge_hash_to_fiber_index.end())
                {
                    auto fiber_ind    = find_iter->second;
                    auto fiber_center = FiberCenters[fiber_ind];
                    auto fiber_radius = FiberRadii[fiber_ind];
                    auto face_conn    = f->get_conn();
                    for (auto i = 0; i < face_conn.size(); i++)
                    {
                        auto coord   = face_conn[i]->get_coord();
                        auto dir_vec = subtract_vec(coord, fiber_center);
                        normalize(dir_vec);
                        coord           = add_vec(fiber_center, scale_vec(dir_vec, fiber_radius));
                        (*face_conn[i]) = coord;
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

        mesh.add_to_self(FibersMesh);
    }

    // If bad meshes are commonly created due to duplicate nodes (you can know by the boundary mesh
    // that is output when you type "Debug" after the automatic periodic constraints command in
    // Beta), then consider making the compare tolerance a function of the shortest edge length.
    MeshManipulation<MeshType>::remove_duplicate_nodes(mesh);
    mesh.update_global_numbering();
}

template<typename MeshType>
inline MeshType MeshCreationFiberMatrix<MeshType>::extract_boundary_mesh(
  std::vector<std::vector<double>>& NodeCoords,
  std::vector<std::vector<int>>& BoundarySegments,
  std::vector<int>& SegmentMarkers)
{
    MeshType tmp_mesh;
    std::vector<NodePtr> nodes;
    for (auto n : NodeCoords)
    {
        auto nPtr = tmp_mesh.add_node_by_coord(n);
        nodes.push_back(nPtr);
    }
    for (auto i = 0; i < BoundarySegments.size(); i++)
    {
        auto b    = BoundarySegments[i];
        auto ePtr = tmp_mesh.add_element();
        ePtr->initialize(ElemGeomType::L1D2N, { nodes[b[0]], nodes[b[1]] });
        if (SegmentMarkers[i] < 0)
            ePtr->set_domain_id(SegmentMarkers[i] + 100);
    }
    return tmp_mesh;
}
*/
}