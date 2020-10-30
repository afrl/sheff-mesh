#include "mesh/external_interfaces/fea_io/TypeConversions.h"

namespace mesh
{
LCSAssignmentType mesh::convert_lcs_assignment(fea_io::LCSAssignment lcs_assignment)
{
    if (lcs_assignment == fea_io::LCSAssignment::Component)
        return LCSAssignmentType::Component;
    else if (lcs_assignment == fea_io::LCSAssignment::Element)
        return LCSAssignmentType::Element;
    else if (lcs_assignment == fea_io::LCSAssignment::Node)
        return LCSAssignmentType::Node;
    else if (lcs_assignment == fea_io::LCSAssignment::None)
        return LCSAssignmentType::None;
    else
        throw _UTIL_MSG_EXCEPT("Cannot convert from fea_io::LCSAssignment.");
}

fea_io::LCSAssignment convert_lcs_assignment(LCSAssignmentType lcs_assignment)
{
    if (lcs_assignment == LCSAssignmentType::Component)
        return fea_io::LCSAssignment::Component;
    else if (lcs_assignment == LCSAssignmentType::Element)
        return fea_io::LCSAssignment::Element;
    else if (lcs_assignment == LCSAssignmentType::Node)
        return fea_io::LCSAssignment::Node;
    else if (lcs_assignment == LCSAssignmentType::None)
        return fea_io::LCSAssignment::None;
    else
        throw _UTIL_MSG_EXCEPT("Cannot convert from fea_io::LCSAssignment.");
}

ElementType convert_element_type(fea_io::ElementType type)
{
    // 1D
    if (type == fea_io::ElementType::L1D2N)
        return ElementType::L1D2N;
    else if (type == fea_io::ElementType::L1D3N)
        return ElementType::L1D3N;
    // 2D
    else if (type == fea_io::ElementType::T2D3N)
        return ElementType::T2D3N;
    else if (type == fea_io::ElementType::Q2D4N)
        return ElementType::Q2D4N;
    else if (type == fea_io::ElementType::C2D4N)
        return ElementType::C2D4N;
    else if (type == fea_io::ElementType::T2D6N)
        return ElementType::T2D6N;
    else if (type == fea_io::ElementType::Q2D8N)
        return ElementType::Q2D8N;
    else if (type == fea_io::ElementType::C2D6N)
        return ElementType::C2D6N;
    // 3D Linear
    else if (type == fea_io::ElementType::T3D4N)
        return ElementType::T3D4N;
    else if (type == fea_io::ElementType::P3D5N)
        return ElementType::P3D5N;
    else if (type == fea_io::ElementType::W3D6N)
        return ElementType::W3D6N;
    else if (type == fea_io::ElementType::H3D8N)
        return ElementType::H3D8N;
    else if (type == fea_io::ElementType::Q3D4N)
        return ElementType::Q3D4N;
    else if (type == fea_io::ElementType::C3D8N)
        return ElementType::C3D8N;
    // 3D Quadratic
    else if (type == fea_io::ElementType::T3D10N)
        return ElementType::T3D10N;
    else if (type == fea_io::ElementType::P3D13N)
        return ElementType::P3D13N;
    else if (type == fea_io::ElementType::W3D15N)
        return ElementType::W3D15N;
    else if (type == fea_io::ElementType::H3D20N)
        return ElementType::H3D20N;
    else if (type == fea_io::ElementType::Q3D8N)
        return ElementType::Q3D8N;
    else if (type == fea_io::ElementType::C3D16N)
        return ElementType::C3D16N;
    throw _UTIL_MSG_EXCEPT("No conversion available.");
}

fea_io::ElementType convert_element_type(ElementType type)
{
    // 1D
    switch (type)
    {
        case ElementType::L1D2N: return fea_io::ElementType::L1D2N;
        case ElementType::L1D3N:
            return fea_io::ElementType::L1D3N;
            // 2D
        case ElementType::T2D3N: return fea_io::ElementType::T2D3N;
        case ElementType::Q2D4N: return fea_io::ElementType::Q2D4N;
        case ElementType::C2D4N: return fea_io::ElementType::C2D4N;
        case ElementType::T2D6N: return fea_io::ElementType::T2D6N;
        case ElementType::Q2D8N: return fea_io::ElementType::Q2D8N;
        case ElementType::C2D6N: return fea_io::ElementType::C2D6N;
        case ElementType::POLYGON:
            return fea_io::ElementType::POLYGON;
            // 3D Linear
        case ElementType::T3D4N: return fea_io::ElementType::T3D4N;
        case ElementType::P3D5N: return fea_io::ElementType::P3D5N;
        case ElementType::W3D6N: return fea_io::ElementType::W3D6N;
        case ElementType::H3D8N: return fea_io::ElementType::H3D8N;
        case ElementType::Q3D4N: return fea_io::ElementType::Q3D4N;
        case ElementType::C3D8N:
            return fea_io::ElementType::C3D8N;
            // 3D Quadratic
        case ElementType::T3D10N: return fea_io::ElementType::T3D10N;
        case ElementType::P3D13N: return fea_io::ElementType::P3D13N;
        case ElementType::W3D15N: return fea_io::ElementType::W3D15N;
        case ElementType::H3D20N: return fea_io::ElementType::H3D20N;
        case ElementType::Q3D8N: return fea_io::ElementType::Q3D8N;
        case ElementType::C3D16N: return fea_io::ElementType::C3D16N;
        case ElementType::POLYHEDRON: return fea_io::ElementType::POLYHEDRON;
        default: throw _UTIL_MSG_EXCEPT("Failed to convert element type to fea_io type.");
    }
}

}
