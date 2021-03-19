#include "mesh/algorithms/LCSAssignments.h"

#include "simple_math/Utility.h"

namespace mesh
{
void assign_lcs(Mesh& m,
                const std::string& element_set_name,
                const simple_math::Point& x_vector,
                const simple_math::Point& y_vector)
{
    auto& group = m.elements().group(element_set_name);
    auto new_q  = simple_math::Quaternion::FromTwoVectors(simple_math::point_as_vector(x_vector),
                                                         simple_math::point_as_vector(y_vector));
    for (auto& e_i : group.element_indices)
    {
        auto e = m.elements()[e_i];
        for (auto& q : e.local_coord_sys())
        {
            q = new_q;
        }
    }
}
}