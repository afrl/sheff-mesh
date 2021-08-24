#pragma once

#include "FaceTypes.h"
#include "mesh/Defines.h"
#include "utility/Exceptions.h"

#include <vector>

namespace mesh
{
/// Dev note: if you add new ones, do not change the number of old types so that old
/// mesh files remain compatible.
/// IMPORTANT NOTE: if you add a new element type that you want to be able to vizualize
/// using ParaView, you need to update the fea_io library, in any related functions.
/// Additionally, you need to update the ParaViewReader project file vtkFEAIOReader.cxx methods
/// get_to_vtk_cell_map and get_to_vtk_reordering.
// TODO check names for ambiguity
enum class ElementType : ConnIndex
{
    UNKNGEOM = 0, ///< Used if a conversion fails
    // 1D types 10-19
    DEF1D = 10, ///< 1-D Default
    L1D2N = 12, ///< 1-D Linear 2-node linear element
    L1D3N = 13, ///< 1-D Linear 3-node linear element
    // 2D types 20-24 is linear and 25-29 is quadratic
    DEF2D   = 20, ///< 2-D Default
    T2D3N   = 21, ///< triangular 2D 3-node element (linear)
    Q2D4N   = 22, ///< quadrilateral 4-node linear element
    C2D4N   = 24, ///< line 4-node cohesive element (linear)
    T2D6N   = 25, ///< triangular 2D 6-node element (quadratic)
    Q2D8N   = 26, ///< quadrilateral 8-node quadratic element (quadratic)
    C2D6N   = 28, ///< line 6-node cohesive element (quadratic)
    POLYGON = 29, /// < general 2D N-sided polygon
    // Linear 3D 30-39
    DEF3D      = 30, ///< 3-D Default
    T3D4N      = 31, ///< 3D tetrahedral 4-node  (linear) element
    P3D5N      = 35, ///< 3D pyramid 5-node  (linear) element
    W3D6N      = 32, ///< 3D wedge 6-node  (linear) element
    H3D8N      = 33, ///< 3D hexahedral 8-node (linear) element
    Q3D4N      = 36, ///< quadrilateral shell element  (linear)
    C3D8N      = 38, ///< quadrilateral 8-node cohesive element (linear)
    POLYHEDRON = 39, ///< general 3D N-sided polyhedron
    // Quadratic 3D 41-49
    T3D10N = 41, ///< 3D tetrahedral 10-node  (quadratic) element
    P3D13N = 45, ///< 3D pyramid 13-node (quadratic) element
    W3D15N = 42, ///< 3D wedge 15-node  (quadratic) element
    H3D20N = 43, ///< 3D hexahedral 20-node (quadratic) element
    Q3D8N  = 46, ///< quadrilateral shell element  (quadratic)
    C3D16N = 48  ///< quadrilateral 16-node cohesive element (quadratic)
};
/// The maximum value that a ElementType enumeration can have.  This is only
/// used for classes like Mapping, which use the underlying type as an index.
static constexpr ConnIndex MAX_ELEMENT_TYPE_VALUE = 48;

enum class ElementTypeCategory : int8_t
{
    SOLID,
    INTERFACIAL,
    ANY ///< Used an option when indicating which category to work with
};

/// Returns the category (solid, interfacial, etc.) of the given element type.
static ElementTypeCategory get_category(const ElementType& type);

static std::string to_string(const ElementType& type);

/// Returns the number of nodes that the geometry type will have.
static ConnIndex get_num_nodes_for_type(const ElementType& geom_type);

/// This method will use the number of dimensions and number of nodes to infer
/// the type of solid element.
static ElementType infer_type(unsigned int dimensionality, ConnIndex num_nodes);

/// For 1D elements, returns 1.  For 2D elements, returns 2. For 3D elements, returns 3.
static unsigned int get_dimensionality(const ElementType& type);

/// Given a geometry type, this function will return the linear version
/// if one exists.  Otherwise it will return UNKNGEOM.
static ElementType get_linearized_type(const ElementType& geom_type);
/// Returns the permuation vector to go from the quadratic connectivity
/// to the linear connectivity.  If the type is already linear, then
/// an empty vector is returned.  Will throw for unrecognized types.
static std::vector<ConnIndex> get_linearized_permutation(const ElementType& geom_type);
/// Returns a bool indicating if the geometry type is linear or not
static bool is_type_linear(const ElementType& geom_type);
/// Returns the element order (1 - linear, 2 - quadratic, etc.)
static int get_order(const ElementType& geom_type);

/// Returns the corresponding quadratic type for a given linear type
static ElementType get_quadratic_type(const ElementType& geom_type);
/// For a given type of linear element, this will return the connectivity index
/// map from linear to quadratic for the nodes that already exist (corner nodes)
/// Example for T2D3N, the vector will be:
/// [0]: 0
/// [1]: 2
/// [2]: 4
/// which means that the first node in the linear connectivity belongs at index 0
/// in the quadratic connectivity, etc.
static std::vector<ConnIndex> get_quadratic_indices_of_linear_nodes(const ElementType& geom_type);
/// Given an element geometry type, this function will fill the vectors:
/// - quadratic_indices: the indicides of hte midside nodes in the quadratic connectivity vector
/// - linear_node_pair_indices: the corresponding pair the midside node lies on in terms of the
/// linear
///   connectivity vector.
/// \attention Will overwrite anything in the vectors!
static void get_quadratic_indices_of_midside_nodes(
  const ElementType& geom_type,
  std::vector<ConnIndex>& quadratic_indices,
  std::vector<std::pair<ConnIndex, ConnIndex>>& linear_node_pair_indices);

/// Given an element geometry type, returns the edges of the element
const std::vector<EdgeInfo> get_edge_info(const ElementType& geom_type);
/// Given an element geometry type, returns the faces of the element
const std::vector<FaceInfo>& get_face_info(const ElementType& geom_type);
/// Given an element geometry type, returns the faces of the element as if each
/// higher-order edge is a collection of linear 2-node edges, and the nodes of
/// each face are ordered clockwise.
static const std::vector<FacetInfo> get_facet_info(const ElementType& geom_type);
/// Returns the parametric coordinates of each of the nodes in the element
static const std::vector<Point> get_node_parametric_coordinates(const ElementType& geom_type);
/// Returns the element type that is equivalent to the given face type.  Throws if the face cannot
/// be converted.
static ElementType element_equivalent_for_face(const FaceType& face_type);

/////////////
// Helpers //
/////////////

/// A helper used in get_edge_info and get_face_info
static std::vector<std::vector<ConnIndex>> connectivity_helper(const ConnIndex num_faces,
                                                               const ConnIndex* nodes_per_face,
                                                               const ConnIndex* face_nodes);

//////////////////
// Definitiions //
//////////////////

namespace detail
{
static const std::vector<FacetInfo> get_facet_info_helper_linear(
  const std::vector<FaceInfo>& face_info)
{
    std::vector<FacetInfo> facet_info;
    facet_info.reserve(face_info.size());
    for (auto& face : face_info)
    {
        facet_info.push_back(
          { face.connectivity, std::vector<int8_t>(face.edges.begin(), face.edges.end()) });
    }
    return facet_info;
}

static const std::vector<FacetInfo> get_facet_info_helper_quadratic(
  const std::vector<FaceInfo>& face_info)
{
    std::vector<FacetInfo> facet_info;
    facet_info.reserve(face_info.size());
    for (auto& face : face_info)
    {
        std::vector<int8_t> segment_markers;
        segment_markers.reserve(2 * face.edges.size());
        for (auto& edge_index : face.edges)
        {
            segment_markers.push_back(static_cast<int8_t>(edge_index));
            segment_markers.push_back(static_cast<int8_t>(edge_index));
        }

        facet_info.push_back({ face.connectivity, segment_markers });
    }
    return facet_info;
}
}

const std::vector<FacetInfo> get_facet_info(const ElementType& geom_type)
{
    auto face_info = get_face_info(geom_type);

    switch (geom_type)
    {
        // For these linear 3D elements, the faces for all of the following types are
        // already ordered clock-wise, so no modification to the connectivity from
        // face_info is needed.  Additionaly, each edge already has 2-nodes, so the
        // edges and segments are equivalent (no modification needed).
        case ElementType::T3D4N:
        case ElementType::P3D5N:
        case ElementType::W3D6N:
        case ElementType::H3D8N: return detail::get_facet_info_helper_linear(face_info);
        // For these quadratic 3D elements, the faces for all of the following types are
        // already ordered clock-wise, so no modification to the connectivity from
        // face_info is needed.  However, each edge has 3-nodes, so the edge index
        // needs to be repeated for each edge.
        case ElementType::T3D10N:
        case ElementType::P3D13N:
        case ElementType::W3D15N:
        case ElementType::H3D20N: return detail::get_facet_info_helper_quadratic(face_info);
        case ElementType::C3D8N:
        case ElementType::C3D16N:
        default: throw _UTIL_MSG_EXCEPT("Facet converion is not implemented.");
    }
}

const std::vector<Point> get_node_parametric_coordinates(const ElementType& geom_type)
{
    switch (geom_type)
    {
        // 1D
        case ElementType::L1D2N: return { { -1.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 } }; break;
        case ElementType::L1D3N:
            return { { -1.0, 0.0, 0.0 }, { 0.5, 0.0, 0.0 }, { 1.0, 0.0, 0.0 } };
            break;
        // 2D
        case ElementType::T2D3N:
            return { { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 } };
            break;
        case ElementType::Q2D4N:
            return { { -1.0, -1.0, 0.0 },
                     { 1.0, -1.0, 0.0 },
                     { 1.0, 1.0, 0.0 },
                     { -1.0, 1.0, 0.0 } };
            break;
        case ElementType::T2D6N:
            return { { 0.0, 0.0, 0.0 }, { 0.5, 0.0, 0.0 }, { 1.0, 0.0, 0.0 },
                     { 0.5, 0.5, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.5, 0.0 } };
            break;
        case ElementType::Q2D8N:
            return {
                { -1.0, -1.0, 0.0 }, { 0.0, -1.0, 0.0 }, { 1.0, -1.0, 0.0 }, { 1.0, 0.0, 0.0 },
                { 1.0, 1.0, 0.0 },   { 0.0, 1.0, 0.0 },  { -1.0, 1.0, 0.0 }, { -1.0, 0.0, 0.0 }
            };
            break;
        // 3D
        case ElementType::T3D4N:
            return { { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };
            break;
        case ElementType::P3D5N:
            return { { -1.0, -1.0, -1.0 },
                     { 1.0, -1.0, -1.0 },
                     { 1.0, 1.0, -1.0 },
                     { -1.0, 1.0, -1.0 },
                     { 0.0, 0.0, 1.0 } };
            break;
        case ElementType::W3D6N:
            return { { 0.0, 0.0, -1.0 }, { 1.0, 0.0, -1.0 }, { 0.0, 1.0, -1.0 },
                     { 0.0, 0.0, 1.0 },  { 1.0, 0.0, 1.0 },  { 0.0, 1.0, 1.0 } };
            break;
        case ElementType::H3D8N:
            return { { -1.0, -1.0, -1.0 }, { -1.0, 1.0, -1.0 }, { -1.0, 1.0, 1.0 },
                     { -1.0, -1.0, 1.0 },  { 1.0, -1.0, -1.0 }, { 1.0, 1.0, -1.0 },
                     { 1.0, 1.0, 1.0 },    { 1.0, -1.0, 1.0 } };
            break;
        case ElementType::T3D10N:
            return { { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 },
                     { 0.5, 0.0, 0.0 }, { 0.5, 0.5, 0.0 }, { 0.0, 0.5, 0.0 }, { 0.0, 0.0, 0.5 },
                     { 0.5, 0.0, 0.5 }, { 0.0, 0.5, 0.5 } };
            break;
        case ElementType::P3D13N:
            return { { -1.0, -1.0, -1.0 }, { 1.0, -1.0, -1.0 }, { 1.0, 1.0, -1.0 },
                     { -1.0, 1.0, -1.0 },  { 0.0, 0.0, 1.0 },   { 0.0, -1.0, -1.0 },
                     { 1.0, 0.0, -1.0 },   { 0.0, 1.0, -1.0 },  { -1.0, 0.0, -1.0 },
                     { 0.0, -1.0, 0.0 },   { 1.0, 0.0, 0.0 },   { 0.0, 1.0, 0.0 },
                     { -1.0, 0.0, 0.0 } };
            break;
        case ElementType::W3D15N:
            return { { 0.0, 0.0, -1.0 }, { 0.5, 0.0, -1.0 }, { 1.0, 0.0, -1.0 }, { 0.5, 0.5, -1.0 },
                     { 0.0, 1.0, -1.0 }, { 0.0, 0.5, -1.0 }, { 0.0, 0.0, 0.0 },  { 1.0, 0.0, 0.0 },
                     { 0.0, 1.0, 0.0 },  { 0.0, 0.0, 1.0 },  { 0.5, 0.0, 1.0 },  { 1.0, 0.0, 1.0 },
                     { 0.5, 0.5, 1.0 },  { 0.0, 1.0, 1.0 },  { 0.0, 0.5, 1.0 } };
            break;
        case ElementType::H3D20N:
            return { { -1.0, -1.0, -1.0 }, { -1.0, 0.0, -1.0 }, { -1.0, 1.0, -1.0 },
                     { -1.0, 1.0, 0.0 },   { -1.0, 1.0, 1.0 },  { -1.0, 0.0, 1.0 },
                     { -1.0, -1.0, 1.0 },  { -1.0, -1.0, 0.0 }, { 0.0, -1.0, -1.0 },
                     { 0.0, 1.0, -1.0 },   { 0.0, 1.0, 1.0 },   { 0.0, -1.0, 1.0 },
                     { 1.0, -1.0, -1.0 },  { 1.0, 0.0, -1.0 },  { 1.0, 1.0, -1.0 },
                     { 1.0, 1.0, 0.0 },    { 1.0, 1.0, 1.0 },   { 1.0, 0.0, 1.0 },
                     { 1.0, -1.0, 1.0 },   { 1.0, -1.0, 0.0 } };
            break;
        default: throw _UTIL_MSG_EXCEPT("Parametric coordinates not implemented for element type.");
    }
}

ElementType element_equivalent_for_face(const FaceType& face_type)
{
    switch (face_type)
    {
        case FaceType::L2D2N: return ElementType::L1D2N;
        case FaceType::L2D3N: return ElementType::L1D3N;
        case FaceType::T3D3N: return ElementType::T2D3N;
        case FaceType::Q3D4N: return ElementType::Q2D4N;
        case FaceType::T3D6N: return ElementType::T2D6N;
        case FaceType::Q3D8N: return ElementType::Q2D8N;
    }
    throw _UTIL_MSG_EXCEPT("There is no element equivalent for the given face type.");
}

std::vector<std::vector<ConnIndex>> connectivity_helper(const ConnIndex num_faces,
                                                        const ConnIndex* nodes_per_face,
                                                        const ConnIndex* face_nodes)
{
    std::vector<std::vector<ConnIndex>> faces_connectivities(num_faces);
    int FaceNodesIndex = 0;
    for (int face = 0; face < num_faces; face++)
    {
        faces_connectivities[face].resize(nodes_per_face[face]);
        for (int node = 0; node < nodes_per_face[face]; node++)
        {
            faces_connectivities[face][node] = face_nodes[FaceNodesIndex];
            FaceNodesIndex++;
        }
    }
    return faces_connectivities;
}

ElementTypeCategory get_category(const ElementType& type)
{
    switch (type)
    {
        case ElementType::DEF1D:
        case ElementType::L1D2N:
        case ElementType::L1D3N: return ElementTypeCategory::SOLID;
        // 2D types 20-24 is linear and 25-29 is quadratic
        case ElementType::DEF2D:
        case ElementType::T2D3N:
        case ElementType::Q2D4N:
        case ElementType::T2D6N:
        case ElementType::Q2D8N: return ElementTypeCategory::SOLID;
        case ElementType::POLYGON:
        case ElementType::C2D4N:
        case ElementType::C2D6N: return ElementTypeCategory::INTERFACIAL;
        // 3D
        case ElementType::DEF3D:
        case ElementType::T3D4N:
        case ElementType::P3D5N:
        case ElementType::W3D6N:
        case ElementType::H3D8N:
        case ElementType::Q3D4N:
        case ElementType::POLYHEDRON:
        case ElementType::T3D10N:
        case ElementType::P3D13N:
        case ElementType::W3D15N:
        case ElementType::H3D20N:
        case ElementType::Q3D8N: return ElementTypeCategory::SOLID;
        case ElementType::C3D8N:
        case ElementType::C3D16N: return ElementTypeCategory::INTERFACIAL;
    }
    throw _UTIL_MSG_EXCEPT("A category is not known for element type: " + to_string(type));
}

std::string to_string(const ElementType& type)
{
    switch (type)
    {
        case ElementType::UNKNGEOM: return "UNKNGEOM";
        case ElementType::DEF1D: return "DEF1D";
        case ElementType::L1D2N: return "L1D2N";
        case ElementType::L1D3N: return "L1D3N";
        // 2D types 20-24 is linear and 25-29 is quadratic
        case ElementType::DEF2D: return "DEF2D";
        case ElementType::T2D3N: return "T2D3N";
        case ElementType::Q2D4N: return "Q2D4N";
        case ElementType::C2D4N: return "C2D4N";
        case ElementType::T2D6N: return "T2D6N";
        case ElementType::Q2D8N: return "Q2D8N";
        case ElementType::C2D6N: return "C2D6N";
        case ElementType::POLYGON: return "POLYGON";
        // Linear 3D 30-39
        case ElementType::DEF3D: return "DEF3D";
        case ElementType::T3D4N: return "T3D4N";
        case ElementType::P3D5N: return "P3D5N";
        case ElementType::W3D6N: return "W3D6N";
        case ElementType::H3D8N: return "H3D8N";
        case ElementType::Q3D4N: return "Q3D4N";
        case ElementType::C3D8N: return "C3D8N";
        case ElementType::POLYHEDRON: return "POLYHEDRON";
        // Quadratic 3D 41-49
        case ElementType::T3D10N: return "T3D10N";
        case ElementType::P3D13N: return "P3D13N";
        case ElementType::W3D15N: return "W3D15N";
        case ElementType::H3D20N: return "H3D20N";
        case ElementType::Q3D8N: return "Q3D8N";
        case ElementType::C3D16N: return "C3D16N";
    }
    throw _UTIL_MSG_EXCEPT("String conversion not implemented for element type.");
}

ConnIndex get_num_nodes_for_type(const ElementType& geom_type)
{
    // 1D
    if (geom_type == ElementType::L1D2N)
        return 2;
    else if (geom_type == ElementType::L1D3N)
        return 3;
    // 2D
    else if (geom_type == ElementType::T2D3N)
        return 3;
    else if (geom_type == ElementType::Q2D4N)
        return 4;
    else if (geom_type == ElementType::C2D4N)
        return 4;
    else if (geom_type == ElementType::T2D6N)
        return 6;
    else if (geom_type == ElementType::Q2D8N)
        return 8;
    else if (geom_type == ElementType::C2D6N)
        return 6;
    // 3D
    else if (geom_type == ElementType::T3D4N)
        return 4;
    else if (geom_type == ElementType::P3D5N)
        return 5;
    else if (geom_type == ElementType::W3D6N)
        return 6;
    else if (geom_type == ElementType::H3D8N)
        return 8;
    else if (geom_type == ElementType::Q3D4N)
        return 4;
    else if (geom_type == ElementType::C3D8N)
        return 8;
    else if (geom_type == ElementType::T3D10N)
        return 10;
    else if (geom_type == ElementType::P3D13N)
        return 13;
    else if (geom_type == ElementType::W3D15N)
        return 15;
    else if (geom_type == ElementType::H3D20N)
        return 20;
    else if (geom_type == ElementType::Q3D8N)
        return 8;
    else if (geom_type == ElementType::C3D16N)
        return 16;
    // No match
    else
        return 0;
}

ElementType infer_type(unsigned int dimensionality, ConnIndex num_nodes)
{

    switch (dimensionality)
    {
        case 1:
            switch (num_nodes)
            {
                case 2: return ElementType::L1D2N;
                case 3: return ElementType::L1D3N;
                default:
                    throw _UTIL_MSG_EXCEPT("Cannot infer element type for 1D with " +
                                           std::to_string(num_nodes) + " nodes.");
            }
        case 2:
            switch (num_nodes)
            {
                case 3: return ElementType::T2D3N;
                case 4: return ElementType::Q2D4N;
                case 6: return ElementType::T2D6N;
                case 8: return ElementType::Q2D8N;
                default:
                    throw _UTIL_MSG_EXCEPT("Cannot infer element type for 2D with " +
                                           std::to_string(num_nodes) + " nodes.");
            }
        case 3:
            switch (num_nodes)
            {
                case 4: return ElementType::T3D4N;
                case 5: return ElementType::P3D5N;
                case 6: return ElementType::W3D6N;
                case 8: return ElementType::H3D8N;
                case 10: return ElementType::T3D10N;
                case 13: return ElementType::P3D13N;
                case 15: return ElementType::W3D15N;
                case 20: return ElementType::H3D20N;
                default:
                    throw _UTIL_MSG_EXCEPT("Cannot infer element type for 3D with " +
                                           std::to_string(num_nodes) + " nodes.");
            }
        default: throw _UTIL_MSG_EXCEPT("Invalid dimensionality.");
    }
}

unsigned int get_dimensionality(const ElementType& type)
{
    using utype = std::underlying_type<ElementType>::type;
    if (static_cast<utype>(type) < 20)
    {
        return 1;
    }
    else if (static_cast<utype>(type) < 30)
    {
        return 2;
    }
    else if (static_cast<utype>(type) < 50)
    {
        return 3;
    }
    throw _UTIL_MSG_EXCEPT("Element type not accounted for in cases.");
}

ElementType get_linearized_type(const ElementType& geom_type)
{
    switch (geom_type)
    {
        case ElementType::L1D3N: return ElementType::L1D2N;
        case ElementType::T2D6N: return ElementType::T2D3N;
        case ElementType::Q2D8N: return ElementType::Q2D4N;
        case ElementType::C2D6N: return ElementType::C2D4N;
        case ElementType::T3D10N: return ElementType::T3D4N;
        case ElementType::P3D13N: return ElementType::P3D5N;
        case ElementType::W3D15N: return ElementType::W3D6N;
        case ElementType::H3D20N: return ElementType::H3D8N;
        case ElementType::Q3D8N: return ElementType::Q3D4N;
        case ElementType::C3D16N: return ElementType::C3D8N;
    }
    return geom_type;
}

bool is_type_linear(const ElementType& geom_type)
{
    // 1D
    if (geom_type == ElementType::L1D2N)
        return true;
    // 2D
    else if (geom_type == ElementType::T2D3N)
        return true;
    else if (geom_type == ElementType::Q2D4N)
        return true;
    else if (geom_type == ElementType::C2D4N)
        return true;
    // 3D
    else if (geom_type == ElementType::T3D4N)
        return true;
    else if (geom_type == ElementType::P3D5N)
        return true;
    else if (geom_type == ElementType::W3D6N)
        return true;
    else if (geom_type == ElementType::H3D8N)
        return true;
    else if (geom_type == ElementType::Q3D4N)
        return true;
    else if (geom_type == ElementType::C3D8N)
        return true;

    // No match
    return false;
}

int get_order(const ElementType& geom_type)
{
    switch (geom_type)
    {
        case ElementType::L1D2N:
        case ElementType::T2D3N:
        case ElementType::Q2D4N:
        case ElementType::T3D4N:
        case ElementType::P3D5N:
        case ElementType::W3D6N:
        case ElementType::H3D8N:
        case ElementType::Q3D4N:
        case ElementType::C3D8N: return 1;
        case ElementType::L1D3N:
        case ElementType::T2D6N:
        case ElementType::Q2D8N:
        case ElementType::C2D6N:
        case ElementType::T3D10N:
        case ElementType::P3D13N:
        case ElementType::W3D15N:
        case ElementType::H3D20N:
        case ElementType::Q3D8N:
        case ElementType::C3D16N: return 2;
        default:
            throw _UTIL_MSG_EXCEPT(
              "Order of element type has not been added to switch or is undefined.");
    }
}

ElementType get_quadratic_type(const ElementType& geom_type)
{
    // 2D cases
    if (geom_type == ElementType::T2D3N)
        return ElementType::T2D6N;
    else if (geom_type == ElementType::Q2D4N)
        return ElementType::Q2D8N;
    else if (geom_type == ElementType::C2D4N)
        return ElementType::C2D6N;
    // 3D cases
    else if (geom_type == ElementType::T3D4N)
        return ElementType::T3D10N;
    else if (geom_type == ElementType::P3D5N)
        return ElementType::P3D13N;
    else if (geom_type == ElementType::W3D6N)
        return ElementType::W3D15N;
    else if (geom_type == ElementType::H3D8N)
        return ElementType::H3D20N;
    else if (geom_type == ElementType::C3D8N)
        return ElementType::C3D16N;
    else
        throw _UTIL_MSG_EXCEPT("Quadratic mapping not implemented for type " +
                               to_string(geom_type));
}

std::vector<ConnIndex> get_quadratic_indices_of_linear_nodes(const ElementType& geom_type)
{
    // 2D cases
    if (geom_type == ElementType::T2D3N)
        return { 0, 2, 4 };
    else if (geom_type == ElementType::Q2D4N)
        return { 0, 2, 4, 6 };
    else if (geom_type == ElementType::C2D4N)
        return { 0, 2, 3, 5 };
    // 3D cases
    else if (geom_type == ElementType::T3D4N)
        return { 0, 1, 2, 3 };
    else if (geom_type == ElementType::P3D5N)
        return { 0, 1, 2, 3, 4 };
    else if (geom_type == ElementType::W3D6N)
        return { 0, 2, 4, 9, 11, 13 };
    else if (geom_type == ElementType::H3D8N)
        return { 0, 2, 4, 6, 12, 14, 16, 18 };
    else if (geom_type == ElementType::C3D8N)
        return { 0, 2, 4, 6, 8, 10, 12, 14 };
    else
        throw _UTIL_MSG_EXCEPT("Quadratic mapping not implemented for type " +
                               to_string(geom_type));
}

void get_quadratic_indices_of_midside_nodes(
  const ElementType& geom_type,
  std::vector<ConnIndex>& quadratic_indices,
  std::vector<std::pair<ConnIndex, ConnIndex>>& linear_node_pair_indices)
{
    // 2D cases
    if (geom_type == ElementType::T2D3N)
    {
        quadratic_indices        = { 1, 3, 5 };
        linear_node_pair_indices = { { 0, 2 }, { 2, 4 }, { 4, 0 } };
    }
    else if (geom_type == ElementType::Q2D4N)
    {
        quadratic_indices        = { 1, 3, 5, 7 };
        linear_node_pair_indices = { { 0, 2 }, { 2, 4 }, { 4, 6 }, { 6, 0 } };
    }
    else if (geom_type == ElementType::C2D4N)
    {
        quadratic_indices        = { 1, 4 };
        linear_node_pair_indices = { { 0, 2 }, { 3, 5 } };
    }
    // 3D cases
    else if (geom_type == ElementType::T3D4N)
    {
        quadratic_indices        = { 4, 5, 6, 7, 8, 9 };
        linear_node_pair_indices = { { 0, 1 }, { 1, 2 }, { 2, 0 }, { 0, 3 }, { 1, 3 }, { 2, 3 } };
    }
    else if (geom_type == ElementType::P3D5N)
    {
        quadratic_indices        = { 5, 6, 7, 8, 9, 10, 11, 12 };
        linear_node_pair_indices = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 },
                                     { 0, 4 }, { 1, 4 }, { 2, 4 }, { 3, 4 } };
    }
    else if (geom_type == ElementType::W3D6N)
    {
        quadratic_indices        = { 1, 3, 5, 6, 7, 8, 10, 12, 14 };
        linear_node_pair_indices = { { 0, 2 },  { 2, 4 },  { 0, 4 },   { 0, 9 }, { 2, 11 },
                                     { 4, 13 }, { 9, 11 }, { 11, 13 }, { 9, 13 } };
    }
    else if (geom_type == ElementType::H3D8N)
    {
        quadratic_indices        = { 1, 3, 5, 7, 8, 9, 10, 11, 13, 15, 17, 19 };
        linear_node_pair_indices = { { 0, 2 },   { 2, 4 },   { 4, 6 },   { 6, 0 },
                                     { 0, 12 },  { 2, 14 },  { 4, 16 },  { 6, 18 },
                                     { 12, 14 }, { 14, 16 }, { 16, 18 }, { 18, 12 } };
    }
    else if (geom_type == ElementType::C3D8N)
    {
        quadratic_indices        = { 1, 3, 5, 7, 9, 11, 13, 15 };
        linear_node_pair_indices = { { 0, 2 },  { 2, 4 },   { 4, 6 },   { 6, 0 },
                                     { 8, 10 }, { 10, 12 }, { 12, 14 }, { 14, 8 } };
    }
    else
    {
        throw _UTIL_MSG_EXCEPT("Quadratic mapping not implemented for type " +
                               to_string(geom_type));
    }
}
};
