#include "mesh/algorithms/UnitElements.h"

#include "mesh/common/QuadratureHelpers.h"
#include "mesh/element/NodeParametricCoordinates.h"
#include "simple_math/Utility.h"

namespace mesh
{
Mesh get_unit_element(const ElementType& type)
{
    auto nodal_parametric_points =
      get_nodal_parametric_points(type, select_gauss_quadrature_scheme(type).type);

    Mesh mesh;
    for (auto i : utility::range(nodal_parametric_points.rows()))
    {
        mesh.nodes().add({ nodal_parametric_points(i, 0),
                           nodal_parametric_points(i, 1),
                           nodal_parametric_points(i, 2) });
    }
    std::vector<LocalIndex> connectivity(nodal_parametric_points.rows());
    std::iota(connectivity.begin(), connectivity.end(), 0);
    mesh.elements().add(type, connectivity);

    return mesh;
}
}