#pragma once

#include "mesh/Defines.h"
#include "mesh/common/QuaternionContainer.h"
#include "mesh/common/ElementTypes.h"

#ifdef USE_FEAIO
#include "fea_io/mesh/MeshEnums.h"
#endif

namespace mesh
{
/// @brief Convergs LCSAssignment from fea_io type to mesh type
LCSAssignmentType convert_lcs_assignment(fea_io::LCSAssignment lcs_assignment);
/// @brief Convergs LCSAssignment from mesh type to fea_io type
fea_io::LCSAssignment convert_lcs_assignment(LCSAssignmentType lcs_assignment);

/// @brief Convergs ElementType from fea_io type to mesh type
ElementType convert_element_type(fea_io::ElementType type);
/// @brief Convergs ElementType from mesh type to fea_io type
fea_io::ElementType convert_element_type(ElementType type);
}