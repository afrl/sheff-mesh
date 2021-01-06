#include "mesh/algorithms/TetgenHelpers.h"

namespace mesh
{
Mesh get_mesh_from_tetgen_input(tetgen::tetgenio& in)
{
    Mesh mesh;

    // Create nodes
    for (int i = 0; i < in.numberofpoints; i++)
    {
        mesh.nodes().add(
          { in.pointlist[i * 3], in.pointlist[i * 3 + 1], in.pointlist[i * 3 + 2] });
    }

    // Create elements
    tetgen::tetgenio::facet* f;
    tetgen::tetgenio::polygon* p;
    for (int i = 0; i < in.numberoffacets; i++)
    {
        f = &in.facetlist[i];
        for (int j = 0; j < f->numberofpolygons; j++)
        {
            p = &f->polygonlist[j];
            std::vector<LocalIndex> element_connectivity(p->numberofvertices);
            for (int j = 0; j < p->numberofvertices; j++)
            {
                element_connectivity[j] = p->vertexlist[j];
            }
            mesh.elements().add(ElementType::POLYGON, element_connectivity);
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
}