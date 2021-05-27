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

    simple_math::Quaternion local_to_global_q;
    {
        simple_math::Vector3D x(x_vector[0], x_vector[1], x_vector[2]);
        x.normalize();
        simple_math::Vector3D y(y_vector[0], y_vector[1], y_vector[2]);
        y.normalize();
        simple_math::Vector3D z = x.cross(y);
        y                       = z.cross(x);
        simple_math::Matrix3D global_to_local;
        for (auto i : utility::range(3))
        {
            global_to_local(0, i) = x[i];
            global_to_local(1, i) = y[i];
            global_to_local(2, i) = z(i);
        }
        local_to_global_q = simple_math::Quaternion(global_to_local.transpose());
    }

    for (auto& i : group.element_indices)
    {
        // m.elements()[i].local_coord_sys().initialize(LCSAssignmentType::Element);
        m.elements()[i].local_coord_sys().initialize(LCSAssignmentType::Node,
                                                     m.elements()[i].connectivity().size());
        for (auto& q : m.elements()[i].local_coord_sys())
        {
            q = local_to_global_q;
        }
    }
}
}