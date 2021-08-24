#include "mesh/common/ElementTypes.h"

namespace mesh
{
////////////////////////////////////////////////////
// Declaration of template, see specializations below
template<ElementType T>
const std::vector<FaceInfo>& get_face_info_kernel();

template<>
const std::vector<FaceInfo>& get_face_info_kernel<ElementType::T3D4N>()
{
    static std::vector<FaceInfo> face_info = { { FaceType::T3D3N, { 1, 3, 0 }, { 0, 2, 1 } },
                                               { FaceType::T3D3N, { 2, 5, 1 }, { 0, 3, 2 } },
                                               { FaceType::T3D3N, { 3, 5, 4 }, { 1, 2, 3 } },
                                               { FaceType::T3D3N, { 0, 4, 2 }, { 0, 1, 3 } } };
    return face_info;
}

template<>
const std::vector<FaceInfo>& get_face_info_kernel<ElementType::P3D5N>()
{
    static std::vector<FaceInfo> face_info = { { FaceType::Q3D4N, { 3, 2, 1, 0 }, { 0, 3, 2, 1 } },
                                               { FaceType::T3D3N, { 0, 5, 4 }, { 0, 1, 4 } },
                                               { FaceType::T3D3N, { 1, 6, 5 }, { 1, 2, 4 } },
                                               { FaceType::T3D3N, { 2, 7, 6 }, { 2, 3, 4 } },
                                               { FaceType::T3D3N, { 3, 4, 7 }, { 3, 0, 4 } } };
    return face_info;
}

template<>
const std::vector<FaceInfo>& get_face_info_kernel<ElementType::W3D6N>()
{
    static std::vector<FaceInfo> face_info = {
        { FaceType::T3D3N, { 2, 1, 0 }, { 0, 2, 1 } },
        { FaceType::T3D3N, { 6, 7, 8 }, { 3, 4, 5 } },
        { FaceType::Q3D4N, { 0, 4, 6, 3 }, { 0, 1, 4, 3 } },
        { FaceType::Q3D4N, { 3, 8, 5, 2 }, { 0, 3, 5, 2 } },
        { FaceType::Q3D4N, { 1, 5, 7, 4 }, { 1, 2, 5, 4 } }
    };
    return face_info;
}

template<>
const std::vector<FaceInfo>& get_face_info_kernel<ElementType::H3D8N>()
{
    static std::vector<FaceInfo> face_info = {
        { FaceType::Q3D4N, { 3, 2, 1, 0 }, { 0, 3, 2, 1 } },
        { FaceType::Q3D4N, { 4, 11, 7, 3 }, { 0, 4, 7, 3 } },
        { FaceType::Q3D4N, { 0, 5, 8, 4 }, { 0, 1, 5, 4 } },
        { FaceType::Q3D4N, { 8, 9, 10, 11 }, { 4, 5, 6, 7 } },
        { FaceType::Q3D4N, { 1, 6, 9, 5 }, { 1, 2, 6, 5 } },
        { FaceType::Q3D4N, { 2, 7, 10, 6 }, { 2, 3, 7, 6 } }
    };
    return face_info;
}

template<>
const std::vector<FaceInfo>& get_face_info_kernel<ElementType::T3D10N>()
{
    static std::vector<FaceInfo> face_info = {
        { FaceType::T3D6N, { 1, 3, 0 }, { 0, 6, 2, 5, 1, 4 } },
        { FaceType::T3D6N, { 2, 5, 1 }, { 0, 7, 3, 9, 2, 6 } },
        { FaceType::T3D6N, { 3, 5, 4 }, { 1, 5, 2, 9, 3, 8 } },
        { FaceType::T3D6N, { 0, 4, 2 }, { 0, 4, 1, 8, 3, 7 } }
    };
    return face_info;
}

template<>
const std::vector<FaceInfo>& get_face_info_kernel<ElementType::P3D13N>()
{
    static std::vector<FaceInfo> face_info = {
        { FaceType::Q3D8N, { 0, 1, 2, 3 }, { 0, 8, 3, 7, 2, 6, 1, 5 } },
        { FaceType::T3D6N, { 0, 5, 4 }, { 0, 5, 1, 10, 4, 9 } },
        { FaceType::T3D6N, { 1, 6, 5 }, { 1, 6, 2, 11, 4, 10 } },
        { FaceType::T3D6N, { 2, 7, 6 }, { 2, 7, 3, 12, 4, 11 } },
        { FaceType::T3D6N, { 3, 4, 7 }, { 3, 8, 0, 9, 4, 12 } }
    };
    return face_info;
}

template<>
const std::vector<FaceInfo>& get_face_info_kernel<ElementType::W3D15N>()
{
    static std::vector<FaceInfo> face_info = {
        { FaceType::T3D6N, { 0, 2, 1 }, { 0, 5, 4, 3, 2, 1 } },
        { FaceType::T3D6N, { 6, 7, 8 }, { 9, 10, 11, 12, 13, 14 } },
        { FaceType::Q3D8N, { 0, 4, 6, 3 }, { 0, 1, 2, 7, 11, 10, 9, 6 } },
        { FaceType::Q3D8N, { 3, 8, 5, 2 }, { 0, 6, 9, 14, 13, 8, 4, 5 } },
        { FaceType::Q3D8N, { 1, 5, 7, 4 }, { 2, 3, 4, 8, 13, 12, 11, 7 } }
    };
    return face_info;
}

template<>
const std::vector<FaceInfo>& get_face_info_kernel<ElementType::H3D20N>()
{
    static std::vector<FaceInfo> face_info = {
        { FaceType::Q3D8N, { 3, 2, 1, 0 }, { 0, 7, 6, 5, 4, 3, 2, 1 } },
        { FaceType::Q3D8N, { 4, 11, 7, 3 }, { 0, 8, 12, 19, 18, 11, 6, 7 } },
        { FaceType::Q3D8N, { 0, 5, 8, 4 }, { 0, 1, 2, 9, 14, 13, 12, 8 } },
        { FaceType::Q3D8N, { 8, 9, 10, 11 }, { 12, 13, 14, 15, 16, 17, 18, 19 } },
        { FaceType::Q3D8N, { 1, 6, 9, 5 }, { 2, 3, 4, 10, 16, 15, 14, 9 } },
        { FaceType::Q3D8N, { 2, 7, 10, 6 }, { 4, 5, 6, 11, 18, 17, 16, 10 } }
    };
    return face_info;
}

// End get_face_info_kernel specializations
////////////////////////////////////////////////////

const std::vector<FaceInfo>& get_face_info(const ElementType& geom_type)
{
    // All the different face definitions
    // Faces defined CCW from outside element and starting at a corner

    // Faces for solid elements
    // 2D Connectivities
    switch (geom_type)
    {
        case ElementType::T2D3N:
        case ElementType::T2D6N:
        case ElementType::Q2D4N:
        case ElementType::Q2D8N:
        case ElementType::C2D4N:
        case ElementType::C2D6N:
            throw _UTIL_MSG_EXCEPT("Faces do not exist for 2D elements, only edges.");
        // 3D Connectivitives
        case ElementType::T3D4N:
            return get_face_info_kernel<ElementType::T3D4N>();
        case ElementType::P3D5N:
            return get_face_info_kernel<ElementType::P3D5N>();
        case ElementType::W3D6N:
            return get_face_info_kernel<ElementType::W3D6N>();
        case ElementType::H3D8N:
            return get_face_info_kernel<ElementType::H3D8N>();
        case ElementType::T3D10N:
            return get_face_info_kernel<ElementType::T3D10N>();
        case ElementType::P3D13N:
            return get_face_info_kernel<ElementType::P3D13N>();
        case ElementType::W3D15N:
            return get_face_info_kernel<ElementType::W3D15N>();
        case ElementType::H3D20N:
            return get_face_info_kernel<ElementType::H3D20N>();
        case ElementType::C3D8N:
        case ElementType::C3D16N:
        default: throw _UTIL_MSG_EXCEPT("Face connectivity not implemented.");
    }
}











const std::vector<EdgeInfo> get_edge_info(const ElementType& geom_type)
{
    std::vector<EdgeInfo> edge_info;
    switch (geom_type)
    {
        // 2D Connectivities
        case ElementType::T2D3N:
            // Linear triangle
            edge_info = { { EdgeType::E2N, { 0, 1 } },
                          { EdgeType::E2N, { 1, 2 } },
                          { EdgeType::E2N, { 2, 0 } } };
            break;
        case ElementType::T2D6N:
            // Quadratic triangle
            edge_info = { { EdgeType::E3N, { 0, 1, 2 } },
                          { EdgeType::E3N, { 2, 3, 4 } },
                          { EdgeType::E3N, { 4, 5, 0 } } };
            break;
        case ElementType::Q2D4N:
            // Linear quadralateral
            edge_info = { { EdgeType::E2N, { 0, 1 } },
                          { EdgeType::E2N, { 1, 2 } },
                          { EdgeType::E2N, { 2, 3 } },
                          { EdgeType::E2N, { 3, 0 } } };
            break;
        case ElementType::Q2D8N:
            // Quadratric quadralateral
            edge_info = { { EdgeType::E3N, { 0, 1, 2 } },
                          { EdgeType::E3N, { 2, 3, 4 } },
                          { EdgeType::E3N, { 4, 5, 6 } },
                          { EdgeType::E3N, { 6, 7, 0 } } };
            break;
        case ElementType::C2D4N:
            // Linear 2D cohesive element
            edge_info = { { EdgeType::E2N, { 0, 1 } }, { EdgeType::E2N, { 2, 3 } } };
            break;
        case ElementType::C2D6N:
            // Quadratic 2D cohesive element
            edge_info = { { EdgeType::E2N, { 0, 1, 2 } }, { EdgeType::E2N, { 5, 4, 3 } } };
            break;
        // 3D Connectivities
        case ElementType::T3D4N:
            // Linear tetrahedral
            edge_info = { { EdgeType::E2N, { 0, 1 } }, { EdgeType::E2N, { 0, 2 } },
                          { EdgeType::E2N, { 0, 3 } }, { EdgeType::E2N, { 1, 2 } },
                          { EdgeType::E2N, { 1, 3 } }, { EdgeType::E2N, { 2, 3 } } };
            break;
        case ElementType::P3D5N:
            // Linear pyramid
            edge_info = { { EdgeType::E2N, { 0, 1 } }, { EdgeType::E2N, { 1, 2 } },
                          { EdgeType::E2N, { 2, 3 } }, { EdgeType::E2N, { 3, 0 } },
                          { EdgeType::E2N, { 0, 4 } }, { EdgeType::E2N, { 1, 4 } },
                          { EdgeType::E2N, { 2, 4 } }, { EdgeType::E2N, { 3, 4 } } };
            break;
        case ElementType::W3D6N:
            // Linear wedge
            edge_info = { { EdgeType::E2N, { 0, 1 } }, { EdgeType::E2N, { 1, 2 } },
                          { EdgeType::E2N, { 2, 0 } }, { EdgeType::E2N, { 0, 3 } },
                          { EdgeType::E2N, { 1, 4 } }, { EdgeType::E2N, { 2, 5 } },
                          { EdgeType::E2N, { 3, 4 } }, { EdgeType::E2N, { 4, 5 } },
                          { EdgeType::E2N, { 5, 3 } } };
            break;
        case ElementType::H3D8N:
            // Linear hexahedral
            edge_info = { { EdgeType::E2N, { 0, 1 } }, { EdgeType::E2N, { 1, 2 } },
                          { EdgeType::E2N, { 2, 3 } }, { EdgeType::E2N, { 3, 0 } },
                          { EdgeType::E2N, { 0, 4 } }, { EdgeType::E2N, { 1, 5 } },
                          { EdgeType::E2N, { 2, 6 } }, { EdgeType::E2N, { 3, 7 } },
                          { EdgeType::E2N, { 4, 5 } }, { EdgeType::E2N, { 5, 6 } },
                          { EdgeType::E2N, { 6, 7 } }, { EdgeType::E2N, { 7, 4 } } };
            break;
        case ElementType::T3D10N:
            // Quadratic tetrahedral
            edge_info = { { EdgeType::E3N, { 0, 4, 1 } }, { EdgeType::E3N, { 0, 6, 2 } },
                          { EdgeType::E3N, { 0, 7, 3 } }, { EdgeType::E3N, { 1, 5, 2 } },
                          { EdgeType::E3N, { 1, 8, 3 } }, { EdgeType::E3N, { 2, 9, 3 } } };
            break;
        case ElementType::P3D13N:
            // Quadratic Pyramid
            edge_info = {
                { EdgeType::E3N, { 0, 5, 1 } },  { EdgeType::E3N, { 1, 6, 2 } },
                { EdgeType::E3N, { 2, 7, 3 } },  { EdgeType::E3N, { 3, 8, 0 } },
                { EdgeType::E3N, { 0, 9, 4 } },  { EdgeType::E3N, { 1, 10, 4 } },
                { EdgeType::E3N, { 2, 11, 4 } }, { EdgeType::E3N, { 3, 12, 4 } },
            };
            break;
        case ElementType::W3D15N:
            // Quadratic wedge
            edge_info = { { EdgeType::E3N, { 0, 1, 2 } },   { EdgeType::E3N, { 2, 3, 4 } },
                          { EdgeType::E3N, { 4, 5, 6 } },   { EdgeType::E3N, { 0, 6, 9 } },
                          { EdgeType::E3N, { 2, 7, 11 } },  { EdgeType::E3N, { 4, 8, 13 } },
                          { EdgeType::E3N, { 9, 10, 11 } }, { EdgeType::E3N, { 11, 12, 13 } },
                          { EdgeType::E3N, { 13, 14, 9 } } };
            break;
        case ElementType::H3D20N:
            // Quadratic hexahedral
            edge_info = { { EdgeType::E3N, { 0, 1, 2 } },    { EdgeType::E3N, { 2, 3, 4 } },
                          { EdgeType::E3N, { 4, 5, 6 } },    { EdgeType::E3N, { 6, 7, 0 } },
                          { EdgeType::E3N, { 0, 8, 12 } },   { EdgeType::E3N, { 2, 9, 14 } },
                          { EdgeType::E3N, { 4, 10, 16 } },  { EdgeType::E3N, { 6, 11, 18 } },
                          { EdgeType::E3N, { 12, 13, 14 } }, { EdgeType::E3N, { 14, 15, 16 } },
                          { EdgeType::E3N, { 16, 17, 18 } }, { EdgeType::E3N, { 18, 19, 12 } } };
            break;
        default: throw _UTIL_MSG_EXCEPT("Edge connectivity not implemented.");
    }
    return edge_info;
}
}