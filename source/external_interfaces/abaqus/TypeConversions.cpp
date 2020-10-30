#include "mesh/external_interfaces/abaqus/TypeConversions.h"

#include "utility/StringUtils.h"

namespace mesh
{
ElementType convert_abaqus_type(const std::string& abaqus_type)
{
    if (utility::istr_is_equal(abaqus_type, "C3D8R") || utility::istr_is_equal(abaqus_type, "C3D8"))
    {
        return ElementType::H3D8N;
    }
    else if (utility::istr_is_equal(abaqus_type, "C3D4"))
    {
        return ElementType::T3D4N;
    }
    else if (utility::istr_is_equal(abaqus_type, "CPS3"))
    {
        return ElementType::T2D3N;
    }
    else if (utility::istr_is_equal(abaqus_type, "T3D2"))
    {
        return ElementType::L1D2N;
    }
    // Note: if you add new types here, make sure the connectivity ordering matches between
    // this library and Abaqus.  if it does not, modify read_abaqus() to permute the connectivity
    // accordingly.
    throw _UTIL_MSG_EXCEPT("Conversion not written for Abaqus type: " + abaqus_type);
}
}