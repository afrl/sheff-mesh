#pragma once

#include "flex_mesh/common/Common.h"
#include "flex_mesh/common/Defines.h"
#include "flex_mesh/common/GeometryTypes.h"
#include "flex_mesh/mesh/Mesh.h"
#include "flex_mesh/node/Node.h"
#include "flex_mesh/element/Element.h"

#include "utility/Range.h"
#include "second_party_mesh/tetgen/tetgen.h"

#include <vector>

namespace flex_mesh
{
    template<typename MeshType>
    static MeshType get_mesh_from_tetgen_input(tetgen::tetgenio& in);

    template<typename MeshType>
    static MeshType get_mesh_from_tetgen_output(tetgen::tetgenio& out);

    /////////////////
    // Definitions //
    /////////////////

    template<typename MeshType>
    MeshType get_mesh_from_tetgen_input(tetgen::tetgenio& in)
    {
        using NodeVector = std::vector<typename MeshType::NodeTypename*>;

        MeshType mesh;

        // Create nodes
        NodeVector nodes;
        nodes.reserve(in.numberofpoints);
        for (unsigned int i = 0; i < in.numberofpoints; i++)
        {
            auto new_node = mesh.add_node_by_coord({
                in.pointlist[i * 3],
                in.pointlist[i * 3 + 1],
                in.pointlist[i * 3 + 2] });
            nodes.push_back(new_node);
        }

        // Create elements
        tetgen::tetgenio::facet* f;
        tetgen::tetgenio::polygon* p;
        for (unsigned int i = 0; i < in.numberoffacets; i++)
        {
            f = &in.facetlist[i];
            for (unsigned int j = 0; j < f->numberofpolygons; j++)
            {
                p = &f->polygonlist[j];
                NodeVector element_connectivity(p->numberofvertices, nullptr);
                for (unsigned int j = 0; j < p->numberofvertices; j++)
                {
                    element_connectivity[j] = nodes[p->vertexlist[j]];
                }
                auto new_element = mesh.elements().add();
                new_element->initialize(ElementType::POLYGON, element_connectivity);
            }
        }

        mesh.update_global_numbering();
        return mesh;
    }

    template<typename MeshType>
    MeshType get_mesh_from_tetgen_output(tetgen::tetgenio& out)
    {
        using NodeVector = std::vector<typename MeshType::NodeTypename*>;

        MeshType mesh;

        // Create nodes
        NodeVector nodes;
        nodes.reserve(out.numberofpoints);
        for (unsigned int i = 0; i < out.numberofpoints; i++)
        {
            auto new_node = mesh.add_node_by_coord({
                out.pointlist[i * 3],
                out.pointlist[i * 3 + 1],
                out.pointlist[i * 3 + 2] });
            nodes.push_back(new_node);
        }

        // Create elements
        NodeVector element_connectivity(4, nullptr);
        for (unsigned int i = 0; i < out.numberoftetrahedra; i++)
        {
            for (unsigned int j = 0; j < 4; j++)
            {
                element_connectivity[j] = nodes[out.tetrahedronlist[i * 4 + j]];
            }
            auto new_element = mesh.elements().add();
            new_element->initialize(ElementType::T3D4N, element_connectivity);
        }
        mesh.update_global_numbering();
        return mesh;
    }
}