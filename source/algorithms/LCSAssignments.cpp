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

    simple_math::Quaternion local_to_global_q = x_y_to_quaterion(x_vector, y_vector);

    for (auto& i : group.element_indices)
    {
        m.elements()[i].local_coord_sys().initialize(LCSAssignmentType::Element);
        for (auto& q : m.elements()[i].local_coord_sys())
        {
            q = local_to_global_q;
        }
    }
}

void assign_lcs(ElementAccessor& element,
                const simple_math::Point& x_vector,
                const simple_math::Point& y_vector)
{
    simple_math::Quaternion local_to_global_q = x_y_to_quaterion(x_vector, y_vector);
    if (element.local_coord_sys().size() == 0)
    {
        element.local_coord_sys().initialize(LCSAssignmentType::Element);
    }
    for (auto& q : element.local_coord_sys())
    {
        q = local_to_global_q;
    }
}

simple_math::Quaternion x_y_to_quaterion(const simple_math::Point& local_x_vector,
                                         const simple_math::Point& local_y_vector)
{
    simple_math::Vector3D x(local_x_vector[0], local_x_vector[1], local_x_vector[2]);
    x.normalize();
    simple_math::Vector3D y(local_y_vector[0], local_y_vector[1], local_y_vector[2]);
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
    return simple_math::Quaternion(global_to_local.transpose());
}

simple_math::Quaternion x_z_to_quaterion(const simple_math::Point& local_x_vector,
                                         const simple_math::Point& local_z_vector)
{
    simple_math::Vector3D x(local_x_vector[0], local_x_vector[1], local_x_vector[2]);
    x.normalize();
    simple_math::Vector3D z(local_z_vector[0], local_z_vector[1], local_z_vector[2]);
    z.normalize();
    simple_math::Vector3D y = z.cross(x);
    z                       = x.cross(y);
    simple_math::Matrix3D global_to_local;
    for (auto i : utility::range(3))
    {
        global_to_local(0, i) = x[i];
        global_to_local(1, i) = y[i];
        global_to_local(2, i) = z(i);
    }
    return simple_math::Quaternion(global_to_local.transpose());
}
}