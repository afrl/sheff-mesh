#pragma once

#include "mesh/basis/ElementBasis.h"


namespace mesh
{
simple_math::Matrix3D calc_jacobian(const ElementConstAccessor element,
                                    const simple_math::Point& parametric_coordinate)
{
    simple_math::GeneralMatrix derivatives;
    BasisFunctions::derivatives(element.type(), parametric_coordinate, derivatives);

    auto n_nodes = element.connectivity().size();
    simple_math::GeneralMatrix node_coordinates;
    node_coordinates.resize(n_nodes, 3);
    for (auto i : utility::range(static_cast<mesh::ConnIndex>(n_nodes)))
    {
        for (auto j : utility::range(3))
        {
            node_coordinates(i, j) = element.node(i).coord()[j];
        }
    }

    simple_math::Matrix3D jacobian = derivatives.transpose() * node_coordinates;
    return jacobian;
}
}