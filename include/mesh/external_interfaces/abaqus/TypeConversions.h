#pragma once

#include "mesh/common/ElementTypes.h"

#include <string>

namespace mesh
{
/// Returns the element type from our library that is equivalent to the given Abaqus one.
/// Note: Currently, the type conversions supported do not require a connectivity permutation.
/// But be aware that this might be needed for other Abaqus types.
ElementType convert_abaqus_type(const std::string& abaqus_type);

/// TODO create conversion for connectivity if more element types are to be supported
}