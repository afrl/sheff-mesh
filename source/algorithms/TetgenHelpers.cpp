#include "mesh/algorithms/TetgenHelpers.h"

#include "mesh/algorithms/IO.h"
#include "mesh/algorithms/Search.h"
#include "simple_math/Utility.h"

#include <deque>

namespace mesh
{

Mesh get_surface_from_tetgenio(tetgen::tetgenio& io)
{
    Mesh mesh;

    // Create nodes
    for (int i = 0; i < io.numberofpoints; i++)
    {
        mesh.nodes().add({ io.pointlist[i * 3], io.pointlist[i * 3 + 1], io.pointlist[i * 3 + 2] });
    }

    // Create elements
    tetgen::tetgenio::facet* f;
    tetgen::tetgenio::polygon* p;
    for (int i = 0; i < io.numberoffacets; i++)
    {
        f           = &io.facetlist[i];
        auto marker = io.facetmarkerlist[i];
        for (int j = 0; j < f->numberofpolygons; j++)
        {
            p = &f->polygonlist[j];
            std::vector<LocalIndex> element_connectivity(p->numberofvertices);
            for (int j = 0; j < p->numberofvertices; j++)
            {
                element_connectivity[j] = p->vertexlist[j];
            }
            mesh.elements().add(ElementType::POLYGON, element_connectivity);
            mesh.elements().back().component_id() = marker;
        }
    }

    return mesh;
}

Mesh get_volume_from_tetgenio(tetgen::tetgenio& io)
{
    Mesh mesh;

    // Create nodes
    for (int i = 0; i < io.numberofpoints; i++)
    {
        mesh.nodes().add({ io.pointlist[i * 3], io.pointlist[i * 3 + 1], io.pointlist[i * 3 + 2] });
    }

    // Create elements
    std::vector<LocalIndex> element_connectivity(4);
    for (int i = 0; i < io.numberoftetrahedra; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            element_connectivity[j] = io.tetrahedronlist[i * 4 + j];
        }
        mesh.elements().add(ElementType::T3D4N, element_connectivity);
    }

    if (io.numberoftetrahedronattributes == 1)
    {
        for (int i = 0; i < io.numberoftetrahedra; i++)
        {
            mesh.elements()[i].component_id() = io.tetrahedronattributelist[i];
        }
    }

    return mesh;
}

Mesh tetrahedralize_surface(const Mesh& surface_mesh, const bool delaunay)
{
    auto in_ptr = surface_to_tetgen_input(surface_mesh);

    // Create tetrahedral mesh structure
    tetgen::tetgenbehavior flags;
    flags.plc    = 1;
    flags.object = tetgen::tetgenbehavior::POLY;

    flags.nobisect      = 1;
    flags.nomergefacet  = 1;
    flags.nomergevertex = 1;
    flags.quiet         = 0;

    if (delaunay)
    {
        flags.quality        = 0;
        flags.steinerleft    = 0;
        flags.smooth_maxiter = 0;
        flags.noimprove      = 1;
        flags.epsilon        = 1e-20;
        // flags.facet_separate_ang_tol = 0.0;
        // flags.collinear_ang_tol      = 180.0;
    }
    else
    {
        flags.quality = 1;
    }

    // Parameters
    /*
    if (MaxElementVol > 0.0)
        flags.maxvolume = MaxElementVol;
    if (MinRatio > 0.0)
        flags.minratio = MinRatio;
    if (MinDihedral > 0.0)
        flags.mindihedral = MinDihedral;
    if (NumSteinerPoints >= 0)
        flags.steinerleft = NumSteinerPoints;*/

    // The external call
    tetgen::tetgenio out;
    tetgen::tetrahedralize(&flags, in_ptr.get(), &out);

    auto tetrahedral_mesh = get_volume_from_tetgenio(out);

    return tetrahedral_mesh;
}

std::unique_ptr<tetgen::tetgenio> surface_to_tetgen_input(const Mesh& surface_mesh)
{

    std::unique_ptr<tetgen::tetgenio> in_ptr = std::make_unique<tetgen::tetgenio>();
    tetgen::tetgenio& in                     = *in_ptr;

    // Indicate indices should be zero-based
    in.firstnumber = 0;

    // Copy the node coordinates into the tetgenio struct
    in.numberofpoints = static_cast<int>(surface_mesh.nodes().size());
    in.pointlist      = new REAL[in.numberofpoints * 3];
    for (auto n : surface_mesh.nodes())
    {
        for (unsigned int j = 0; j < 3; j++)
        {
            in.pointlist[3 * n.local_index() + j] = n.point(j);
        }
    }

    // Declare a facet for each surface element
    in.numberoffacets  = surface_mesh.elements().size();
    in.facetlist       = new tetgen::tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];

    // Create the facets (each will have one polygon)
    tetgen::tetgenio::facet* f;
    tetgen::tetgenio::polygon* p;
    for (auto& e : surface_mesh.elements())
    {
        auto& i = e.local_index();

        in.facetmarkerlist[i] = e.local_index() + 1;

        auto conn           = e.connectivity();
        f                   = &in.facetlist[i];
        f->numberofpolygons = 1;
        f->polygonlist      = new tetgen::tetgenio::polygon[f->numberofpolygons];
        f->numberofholes    = 0;
        f->holelist         = NULL;
        p                   = &f->polygonlist[0];
        p->numberofvertices = conn.size();
        p->vertexlist       = new int[p->numberofvertices];

        for (auto j : utility::range(conn.size()))
        {
            p->vertexlist[j] = conn[j];
        }
    }

    return std::move(in_ptr);
}

std::unique_ptr<tetgen::tetgenio> volume_to_tetgen_input(const Mesh& volume_mesh)
{
    std::unique_ptr<tetgen::tetgenio> in_ptr = std::make_unique<tetgen::tetgenio>();
    tetgen::tetgenio& in                     = *in_ptr;

    // Indicate indices should be zero-based
    in.firstnumber = 0;

    // Copy the node coordinates into the tetgenio struct
    in.numberofpoints = static_cast<int>(volume_mesh.nodes().size());
    in.pointlist      = new REAL[in.numberofpoints * 3];
    for (auto n : volume_mesh.nodes())
    {
        for (unsigned int j = 0; j < 3; j++)
        {
            in.pointlist[3 * n.local_index() + j] = n.point(j);
        }
    }

    // Copy tetrahedra
    in.numberoftetrahedra            = volume_mesh.elements().size();
    in.tetrahedronlist               = new int[4 * in.numberoftetrahedra];
    in.numberoftetrahedronattributes = 1;
    in.tetrahedronattributelist      = new double[in.numberoftetrahedra];
    for (auto& e : volume_mesh.elements())
    {
        auto& i = e.local_index();

        for (unsigned int j = 0; j < 4; j++)
        {
            in.tetrahedronlist[4 * i + j] = e.connectivity()[j];
        }
        in.tetrahedronattributelist[i] = e.component_id();
    }

    return std::move(in_ptr);
}

Mesh improve_tetrahedralization(const Mesh& volume_mesh, const ImproveTetParameters& params)
{
    auto in_ptr = volume_to_tetgen_input(volume_mesh);

    return improve_tetrahedralization(*in_ptr, params);
}

Mesh improve_tetrahedralization(tetgen::tetgenio& input, const ImproveTetParameters& params)
{

    // Create tetrahedral mesh structure
    tetgen::tetgenbehavior flags;
    flags.object = tetgen::tetgenbehavior::MESH;
    // I think we always want facets and vertices to be merged? We can add it to the parameter
    // struct if it needs to be an option.
    //
    // flags.nomergefacet  = 1;
    // flags.nomergevertex = 1;
    flags.refine   = 1; // reconstructs tetrahedralization and refines
    flags.quiet    = 0;
    flags.quality  = 1;
    flags.diagnose = 1;
    flags.docheck  = 1;
    flags.verbose  = 1;

    if (params.rel_target_edge_length > 0.0)
    {
        auto extent = mesh::find_extent(
          Eigen::Map<simple_math::DynamicDoubleMatrix>(input.pointlist, input.numberofpoints, 3));
        auto bounding_box_diagonal = simple_math::distance(extent.first, extent.second);
        auto target_edge_length    = params.rel_target_edge_length * bounding_box_diagonal;
        // Based on volume of a "regular" tetrahedron, see
        // https://mathworld.wolfram.com/RegularTetrahedron.html
        flags.maxvolume   = sqrt(2.0) / 12.0 * pow(target_edge_length, 3);
        flags.fixedvolume = 1;
    }
    else
    {
        flags.maxvolume   = -1.0;
        flags.fixedvolume = 0;
    }
    flags.minratio    = params.max_radius_edge_ratio;
    flags.mindihedral = params.min_dihedral_angle;
    flags.steinerleft = params.max_inserted_points;

    // The external call
    tetgen::tetgenio out;
    tetgen::tetrahedralize(&flags, &input, &out);

    auto tetrahedral_mesh = mesh::get_volume_from_tetgenio(out);

    return tetrahedral_mesh;
}

double tetrahedron_signed_volume(const Eigen::RowVector3d& p1,
                                 const Eigen::RowVector3d& p2,
                                 const Eigen::RowVector3d& p3,
                                 const Eigen::RowVector3d& p4)
{
    return (p2 - p1).dot((p3 - p1).cross(p4 - p1)) / 6.0;
}

SegmentedTetrahedralization segmented_tetrahedralization(const Mesh& surface_mesh)
{
    auto in_ptr = surface_to_tetgen_input(surface_mesh);

    // Create tetrahedral mesh structure
    tetgen::tetgenbehavior flags;
    flags.plc                    = 1;
    flags.object                 = tetgen::tetgenbehavior::POLY;
    flags.nobisect               = 1;
    flags.nomergefacet           = 1;
    flags.nomergevertex          = 1;
    flags.quiet                  = 0;
    flags.quality                = 0;
    flags.smooth_maxiter         = 0;
    flags.noimprove              = 1;
    flags.neighout               = 2; // > 1 needed to output tet2face map
    flags.facesout               = 1;
    flags.steinerleft            = 0;
    flags.epsilon                = 0.0;
    flags.facet_separate_ang_tol = 0.0;
    flags.collinear_ang_tol      = 180.0;

    // The external call
    tetgen::tetgenio out;
    tetgen::tetrahedralize(&flags, in_ptr.get(), &out);


    // Region segmentation

    std::vector<int> tet_region_ids(out.numberoftetrahedra, -1);
    int next_start = 0;

    std::deque<int> tets_to_search;

    int curr_region = 0;
    int curr_tet    = next_start;
    while (next_start != out.numberoftetrahedra)
    {
        tet_region_ids[curr_tet] = curr_region;
        for (auto i : utility::range(4))
        {
            const auto& neighbor_shared_face = out.tet2facelist[4 * curr_tet + i];
            // Do not allow adjacency through original surface mesh faces
            if (out.trifacemarkerlist[neighbor_shared_face] == 0)
            {
                const auto& neighbor_tet = out.neighborlist[4 * curr_tet + i];
                assert(tet_region_ids[neighbor_tet] == -1 ||
                       tet_region_ids[neighbor_tet] == curr_region);
                if (tet_region_ids[neighbor_tet] == -1)
                {
                    tets_to_search.push_back(neighbor_tet);
                }
            }
        }

        // Find to next tet
        while (tet_region_ids[curr_tet] == curr_region)
        {
            // No tets more for this region, so start a new region search
            if (tets_to_search.empty())
            {
                for (auto i = next_start; i <= out.numberoftetrahedra; i++)
                {
                    next_start = i;
                    if (i < out.numberoftetrahedra && tet_region_ids[i] == -1)
                    {
                        curr_tet = i;
                        break;
                    }
                }
                curr_region++;
                break;
            }
            else
            { // There more tets to tranverse in this region
                curr_tet = tets_to_search.back();
                tets_to_search.pop_back();
            }
        }
    }

    // End region segmentation

    unsigned int n_regions = curr_region;
    std::vector<std::vector<mesh::LocalIndex>> region_to_tet_map(n_regions);
    simple_math::GeneralMatrix region_interior_points =
      simple_math::GeneralMatrix::Zero(n_regions, 3);
    std::vector<double> max_volume_tet_by_region(n_regions, 0.0);

    // Assign the region ID to the component ID
    auto tetrahedral_mesh = get_volume_from_tetgenio(out);
    for (auto e : tetrahedral_mesh.elements())
    {
        e.component_id() = tet_region_ids[e.local_index()];
        region_to_tet_map[e.component_id()].push_back(e.local_index());
        auto volume = tetrahedron_signed_volume(e.node(0).coord(),
                                                e.node(1).coord(),
                                                e.node(2).coord(),
                                                e.node(3).coord());
        if (volume > max_volume_tet_by_region[e.component_id()])
        {
            auto p                                       = mesh::calc_mean_coordinate(e);
            region_interior_points.row(e.component_id()) = simple_math::point_as_vector(p);
            max_volume_tet_by_region[e.component_id()]   = volume;
        }
    }


    // Debug extract surface

    // Find the vertices we needed
    std::vector<int> out_vertex_to_surface_vertex_map(out.numberofpoints, -1);
    std::vector<int> surface_vertex_to_out_vertex_map(out.numberofpoints, -1);
    int next_surface_index = 0;
    for (auto i : utility::range(out.numberoftrifaces))
    {
        if (out.trifacemarkerlist[i] > 0)
        {
            for (auto j : utility::range(3))
            {
                const auto& out_index = out.trifacelist[3 * i + j];
                if (out_vertex_to_surface_vertex_map[out_index] == -1)
                {
                    out_vertex_to_surface_vertex_map[out_index]          = next_surface_index;
                    surface_vertex_to_out_vertex_map[next_surface_index] = out_index;
                    next_surface_index++;
                }
            }
        }
    }
    Mesh tri_surface;
    for (auto i : utility::range(next_surface_index))
    {
        auto out_index = surface_vertex_to_out_vertex_map[i];
        tri_surface.nodes().add({ out.pointlist[3 * out_index],
                                  out.pointlist[3 * out_index + 1],
                                  out.pointlist[3 * out_index + 2] });
    }
    std::vector<mesh::LocalIndex> conn;
    for (auto i : utility::range(out.numberoftrifaces))
    {
        if (out.trifacemarkerlist[i] > 0)
        {
            for (auto j : utility::range(3))
            {
                conn[j] = out_vertex_to_surface_vertex_map[out.trifacelist[3 * i + j]];
            }
            tri_surface.elements().add(ElementType::T2D3N, conn);
        }
    }

    return { tri_surface, tetrahedral_mesh, n_regions, region_to_tet_map, region_interior_points };
}
}