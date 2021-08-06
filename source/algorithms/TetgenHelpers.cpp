#include "mesh/algorithms/TetgenHelpers.h"

#include "mesh/algorithms/IO.h"

namespace mesh
{
Mesh get_mesh_from_tetgen_input(tetgen::tetgenio& in)
{
    Mesh mesh;

    // Create nodes
    for (int i = 0; i < in.numberofpoints; i++)
    {
        mesh.nodes().add({ in.pointlist[i * 3], in.pointlist[i * 3 + 1], in.pointlist[i * 3 + 2] });
    }

    // Create elements
    tetgen::tetgenio::facet* f;
    tetgen::tetgenio::polygon* p;
    for (int i = 0; i < in.numberoffacets; i++)
    {
        f = &in.facetlist[i];
        auto marker = in.facetmarkerlist[i];
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

Mesh get_mesh_from_tetgen_output(tetgen::tetgenio& out)
{
    Mesh mesh;

    // Create nodes
    for (int i = 0; i < out.numberofpoints; i++)
    {
        mesh.nodes().add(
          { out.pointlist[i * 3], out.pointlist[i * 3 + 1], out.pointlist[i * 3 + 2] });
    }

    // Create elements
    std::vector<LocalIndex> element_connectivity(4);
    for (int i = 0; i < out.numberoftetrahedra; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            element_connectivity[j] = out.tetrahedronlist[i * 4 + j];
        }
        mesh.elements().add(ElementType::T3D4N, element_connectivity);
    }
    return mesh;
}

Mesh tetrahedralize_surface(const Mesh& surface_mesh)
{
    tetgen::tetgenio in;
    // Indicate indices should be zero-based
    in.firstnumber                            = 0;
    std::unique_ptr<tetgen::tetgenio> out_ptr = std::make_unique<tetgen::tetgenio>();
    tetgen::tetgenio& out                     = *out_ptr;

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
    for (auto e : surface_mesh.elements())
    {
        auto& i = e.local_index();

        in.facetmarkerlist[i] = 0;

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

    mesh::write_vtk(get_mesh_from_tetgen_input(in), "tetgen_input.vtk", false);

    // Create tetrahedral mesh structure
    tetgen::tetgenbehavior flags;
    // Switches
    flags.plc           = 1;
    flags.object        = tetgen::tetgenbehavior::POLY;
    flags.nobisect      = 1;
    flags.quality       = 1;
    flags.nomergefacet  = 1;
    flags.nomergevertex = 1;
    flags.quiet         = 0;

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
    tetgen::tetrahedralize(&flags, &in, &out);

    auto tetrahedral_mesh = get_mesh_from_tetgen_output(out);

    return tetrahedral_mesh;
}
}