#include "mesh/algorithms/UnitElements.h"
#include "mesh/basis/BasisFunctions.h"
#include "mesh/element/NodeParametricCoordinates.h"

using namespace mesh;

int main()
{
    {
        Mesh m                       = get_unit_P3D5N();
        auto nodal_parametric_points =
          get_nodal_parametric_points(ElementType::P3D5N, simple_math::QuadratureType::Gauss3D);
        std::cout << "nodal_parametric_points = \n" << nodal_parametric_points << std::endl;

        simple_math::GeneralVector basis_values;
        for (auto i : utility::range(m.nodes().size()))
        {
            std::cout << "node " << i << " = \n" << m.nodes()[i].coord() << std::endl;
            for (auto j : utility::range(3))
            {
                
                assert(nodal_parametric_points(i, j) == m.nodes()[i].point(j));
            }

            BasisFunctions::values(m.elements()[0].type(), m.nodes()[i].point(), basis_values);
            for (auto j : utility::range(m.nodes().size()))
            {
                if (i == j) {
                    assert(basis_values(j) == 1.0);
                }
                else
                {
                    assert(basis_values(j) == 0.0);
                }
            }
        }
    }

    // TODO add check for every type of element.

    return 0;
}